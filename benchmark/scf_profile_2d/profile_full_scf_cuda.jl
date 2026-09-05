#!/usr/bin/env julia

"""Profile two production-like MPO Hartree--Fock iterations on CUDA.

The detailed mode synchronizes around every measured GPU phase. It is a
diagnostic benchmark and is intentionally not the timing path used by normal
production calculations.
"""

using CUDA
using Dates
using ITensors
using MPO_MeanField
using TOML

length(ARGS) == 5 || error(
    "usage: profile_full_scf_cuda.jl CAMPAIGN TASK OUTPUT MODE SCF_ITERATIONS",
)
campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
output = abspath(ARGS[3])
mode = Symbol(ARGS[4])
scf_iterations = parse(Int, ARGS[5])
mode in (:baseline, :lean) || error("MODE must be baseline or lean")
scf_iterations > 0 || error("SCF_ITERATIONS must be positive")
ispath(output) && error("refusing to overwrite existing output: $output")

CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
1 <= task_index <= length(campaign.cases) || error("task index is outside campaign")
case = campaign.cases[task_index]
isnothing(case.spectral_bounds) && error("campaign case requires spectral bounds")
base_params = legacy_parameters(case.model, case.representation, case.solver)
base_params isa ParametersSquare || error("profiler requires a square MPO case")
parameter_names = fieldnames(typeof(base_params))
parameter_values = NamedTuple{parameter_names}(
    Tuple(getfield(base_params, name) for name in parameter_names),
)
params = ParametersSquare(; merge(
    parameter_values,
    (scf_max_iterations=scf_iterations,),
)...)

to_device(value) = ITensors.adapt(CUDA.CuArray, value)
to_host(value) = ITensors.cpu(value)
synchronize() = CUDA.synchronize()

csv(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_row(io, values) = (println(io, join(csv.(values), ',')); flush(io))

mkpath(output)
static_io = open(joinpath(output, "static_phases.csv"), "w")
detail_io = open(joinpath(output, "detail_phases.csv"), "w")
iteration_io = open(joinpath(output, "scf_iterations.csv"), "w")
write_row(static_io, ("phase", "elapsed_time_s", "allocations_bytes", "gc_time_s"))
write_row(detail_io, (
    "scf_iteration", "inner_iteration", "phase", "elapsed_time_s",
    "allocations_bytes", "gc_time_s",
))
write_row(iteration_io, (
    "iteration", "initialization_time_s", "purification_time_s",
    "density_to_host_time_s", "mean_field_time_s", "fields_to_device_time_s",
    "device_diagnostics_time_s", "residuals_time_s", "energy_time_s",
    "mixing_time_s", "measured_iteration_time_s", "purification_iterations",
    "rho_chi", "hartree_chi", "fock_chi", "effective_hamiltonian_chi",
    "mixing_method",
))

static_callback = event -> write_row(static_io, (
    event.phase,
    event.elapsed_time_s,
    hasproperty(event, :allocations_bytes) ? event.allocations_bytes : 0,
    hasproperty(event, :gc_time_s) ? event.gc_time_s : 0.0,
))
detail_callback = event -> write_row(detail_io, (
    event.scf_iteration,
    hasproperty(event, :iteration) ? event.iteration : "",
    event.phase,
    event.elapsed_time_s,
    event.allocations_bytes,
    event.gc_time_s,
))
iteration_callback = event -> write_row(iteration_io, (
    event.iteration,
    event.initialization_time_s,
    event.purification_time_s,
    event.density_to_host_time_s,
    event.mean_field_time_s,
    event.fields_to_device_time_s,
    event.device_diagnostics_time_s,
    event.residuals_time_s,
    event.energy_time_s,
    event.mixing_time_s,
    event.measured_iteration_time_s,
    event.purification_iterations,
    event.rho_bond_dimension,
    event.hartree_bond_dimension,
    event.fock_bond_dimension,
    event.effective_hamiltonian_bond_dimension,
    event.mixing_method,
))

free_before = CUDA.free_memory()
started_at = Dates.now()
println("Constructing static operators for mode=$mode")
flush(stdout)
sys = System(
    params;
    static_to_backend=to_device,
    static_phase_callback=static_callback,
    static_phase_synchronize=synchronize,
)
close(static_io)

measure_stationarity = mode == :baseline
record_energy = mode == :baseline
println(
    "Profiling mode=$mode scf_iterations=$scf_iterations " *
    "stationarity=$measure_stationarity energy=$record_energy",
)
flush(stdout)
wall_started = time_ns()
converged = run_scf!(
    sys,
    case.spectral_bounds...;
    verbose=:all,
    purification_method=case.solver.purification,
    square_fock_method=case.solver.square_fock_method,
    sp2_idempotency_tolerance=case.solver.sp2_idempotency_tolerance,
    sp2_trace_tolerance=(
        case.solver.sp2_relative_trace_tolerance * round(Int, params.density * 2^params.L)
    ),
    record_energy=record_energy,
    stable_iterations=max(case.solver.stable_iterations, scf_iterations + 1),
    require_stationarity=false,
    measure_stationarity=measure_stationarity,
    detect_two_cycles=case.solver.detect_two_cycles,
    mixing_method=case.solver.mixing_method,
    pulay_history=case.solver.pulay_history,
    pulay_warmup=case.solver.pulay_warmup,
    pulay_regularization=case.solver.pulay_regularization,
    pulay_coefficient_limit=case.solver.pulay_coefficient_limit,
    pulay_step_limit=case.solver.pulay_step_limit,
    to_gpu=to_device,
    to_cpu=to_host,
    phase_callback=iteration_callback,
    detail_phase_callback=detail_callback,
    phase_synchronize=synchronize,
    purification_cleanup=CUDA.reclaim,
)
synchronize()
wall_time_s = (time_ns() - wall_started) / 1e9
close(detail_io)
close(iteration_io)

diagnostics = scf_diagnostics(sys)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(started_at),
        "finished_at" => string(Dates.now()),
        "campaign" => campaign_file,
        "label" => case.label,
        "mode" => string(mode),
        "backend" => "cuda",
        "device_name" => CUDA.name(CUDA.device()),
        "matrix_dimension" => 2^params.L,
        "maxdim" => params.itensors_maxdim,
        "cutoff" => params.itensors_tol,
        "requested_scf_iterations" => scf_iterations,
        "completed_scf_iterations" => length(diagnostics.history),
        "scf_converged" => converged,
        "scf_termination_reason" => string(diagnostics.termination_reason),
        "measure_stationarity" => measure_stationarity,
        "record_energy" => record_energy,
        "wall_time_s" => wall_time_s,
        "device_free_memory_before_bytes" => free_before,
        "device_free_memory_after_bytes" => CUDA.free_memory(),
    ))
end

println("profile complete: $output")
println("wall_time_s=$wall_time_s")
flush(stdout)
