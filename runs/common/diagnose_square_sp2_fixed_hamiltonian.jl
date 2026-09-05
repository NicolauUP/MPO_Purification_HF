#!/usr/bin/env julia

"""Diagnose a square campaign's first, fixed-Hamiltonian SP2 purification.

This constructs exactly the initial Hamiltonian used by SCF iteration one,
H0 + S, then runs the production SP2 recurrence without entering SCF or
extracting Hartree/Fock fields. Each `iterations.csv` row is one actual SP2
step before its selected polynomial is applied.

Set `MPO_FIXED_SP2_CHECKPOINT` to an unused path to also write the final
host-side MPO state. Such a checkpoint represents a fixed initial-field SP2
projector, *not* a self-consistent HF solution.
"""

using Dates
using ITensors
using ITensorMPS
using MPO_MeanField
using SHA
using TOML

length(ARGS) in (4, 5, 6, 7) || error("usage: diagnose_square_sp2_fixed_hamiltonian.jl CAMPAIGN_FILE TASK_INDEX MAXDIM OUTPUT_DIRECTORY [BACKEND] [SPECTRAL_LOWER SPECTRAL_UPPER]")
campaign_file = abspath(ARGS[1])
task_index = tryparse(Int, ARGS[2])
maxdim = tryparse(Int, ARGS[3])
output = abspath(ARGS[4])
backend = length(ARGS) in (5, 7) ? Symbol(lowercase(ARGS[5])) : :cpu
isnothing(task_index) && error("TASK_INDEX must be an integer")
isnothing(maxdim) && error("MAXDIM must be an integer")
maxdim > 0 || error("MAXDIM must be positive")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside the campaign")
spec = campaign.runs[task_index]
spec.params isa ParametersSquare || error("fixed square SP2 diagnostic requires ParametersSquare")
ispath(output) && error("refusing to overwrite existing output directory: $output")

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
git_revision(root) = try readchomp(`git -C $root rev-parse HEAD`) catch; "unavailable" end

function with_maxdim(params::ParametersSquare, maxdim::Int)
    ParametersSquare(
        L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
        tci_tol=params.tci_tol, itensors_tol=params.itensors_tol,
        itensors_maxdim=maxdim, density=params.density,
        purification_steps=params.purification_steps,
        scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
        scf_max_iterations=params.scf_max_iterations,
    )
end

function mean_bond_dimension(mpo::MPO)
    length(mpo) <= 1 && return 1.0
    dimensions = [dim(commonind(mpo[index], mpo[index + 1])) for index in 1:(length(mpo) - 1)]
    isempty(dimensions) ? 1.0 : sum(dimensions) / length(dimensions)
end

params = with_maxdim(spec.params, maxdim)
campaign_bounds = (Float64(spec.spectral_bounds[1]), Float64(spec.spectral_bounds[2]))
bounds = if length(ARGS) in (6, 7)
    first_bound_argument = length(ARGS) == 7 ? 6 : 5
    lower = tryparse(Float64, ARGS[first_bound_argument])
    upper = tryparse(Float64, ARGS[first_bound_argument + 1])
    isnothing(lower) && error("SPECTRAL_LOWER must be a float")
    isnothing(upper) && error("SPECTRAL_UPPER must be a float")
    isfinite(lower) && isfinite(upper) && lower < upper || error("spectral bounds must be finite and strictly ordered")
    (lower, upper)
else
    campaign_bounds
end

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end
to_device = backend == :cuda ? (value -> ITensors.adapt(CUDA.CuArray, value)) : identity
synchronize_backend() = backend == :cuda ? CUDA.synchronize() : nothing
device_name = backend == :cuda ? CUDA.name(CUDA.device()) : "CPU"
device_total_memory = backend == :cuda ? CUDA.total_memory() : 0
device_free_memory_before = backend == :cuda ? CUDA.free_memory() : 0

function backend_timed(f)
    synchronize_backend()
    measurement = @timed begin
        value = f()
        synchronize_backend()
        value
    end
    return measurement
end
N = 2 ^ params.L
Ne = round(Int, N * params.density)
repo_root = abspath(joinpath(@__DIR__, "..", ".."))
started_at = now(UTC)
mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))

sys = System(params)
if backend == :cuda
    sys.H0 = to_device(sys.H0)
    sys.VH = to_device(sys.VH)
    sys.VF = to_device(sys.VF)
end
initialization = backend_timed() do
    construct_rho_0(
    sys, params, bounds[1], bounds[2]; method=:sp2,
    verify_spectral_bounds=false, to_gpu=to_device,
    )
end
rho = initialization.value
trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)
hermiticity_tolerance = MPO_MeanField._sp2_hermiticity_tolerance(params)
idempotency_tolerance = 1e-3
previous_idempotency = Inf
stagnant_steps = 0
converged = false
termination_reason = :max_iterations
last_trace = NaN
last_trace_error = NaN
last_idempotency = NaN
last_hermiticity = NaN
last_iteration = 0

open(joinpath(output, "iterations.csv"), "w") do io
    write_csv_row(io, (
        "iteration", "branch", "trace", "trace_squared", "trace_hole",
        "trace_error", "idempotency_residual", "hermiticity_residual",
        "rho_max_chi", "rho_squared_max_chi", "rho_mean_chi",
        "rho_squared_mean_chi", "cap_reached", "step_time_s",
        "step_allocations_bytes", "step_gc_time_s",
    ))
    for iteration in 1:params.purification_steps
        global rho, previous_idempotency, stagnant_steps, converged
        global termination_reason, last_trace, last_trace_error
        global last_idempotency, last_hermiticity, last_iteration
        step = backend_timed() do
            rho_squared = apply(rho, rho;
                cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
            )
            trace_value = real(tr(rho))
            trace_squared = real(tr(rho_squared))
            trace_hole = 2trace_value - trace_squared
            idempotency = MPO_MeanField.idempotency_residual(rho, rho_squared)
            hermiticity = MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
            branch = MPO_MeanField._sp2_branch(trace_value, trace_squared, Ne, trace_tolerance)
            (rho_squared, trace_value, trace_squared, trace_hole, idempotency, hermiticity, branch)
        end
        rho_squared, trace_value, trace_squared, trace_hole, idempotency, hermiticity, branch = step.value
        rho_chi = maxlinkdim(rho)
        rho_squared_chi = maxlinkdim(rho_squared)
        trace_error = abs(trace_value - Ne)
        write_csv_row(io, (
            iteration, branch, trace_value, trace_squared, trace_hole, trace_error,
            idempotency, hermiticity, rho_chi, rho_squared_chi,
            mean_bond_dimension(rho), mean_bond_dimension(rho_squared),
            rho_chi >= maxdim || rho_squared_chi >= maxdim,
            step.time, step.bytes, step.gctime,
        ))
        flush(io)

        last_trace = trace_value
        last_trace_error = trace_error
        last_idempotency = idempotency
        last_hermiticity = hermiticity
        last_iteration = iteration
        if trace_error <= trace_tolerance && idempotency < idempotency_tolerance && hermiticity <= hermiticity_tolerance
            converged = true
            termination_reason = :idempotency_threshold
            break
        end
        if idempotency >= previous_idempotency * (1 - 1e-12)
            stagnant_steps += 1
        else
            stagnant_steps = 0
        end
        if stagnant_steps >= 8
            termination_reason = :stagnation
            break
        end
        previous_idempotency = idempotency
        rho = branch == :square ? rho_squared : +(
            2.0 * rho, -rho_squared;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
    end
end

device_free_memory_after = backend == :cuda ? CUDA.free_memory() : 0
if backend == :cuda
    open(joinpath(output, "cuda_status.txt"), "w") do io
        CUDA.versioninfo(io)
        println(io)
        CUDA.pool_status(io)
    end
end

# A checkpoint is opt-in so the established diagnostic remains lightweight.
# HDF5 serialization is host-side; move only the final compressed state back
# after the timed purification has completed.
checkpoint_path = get(ENV, "MPO_FIXED_SP2_CHECKPOINT", "")
if !isempty(checkpoint_path)
    ispath(checkpoint_path) && error("refusing to overwrite checkpoint: $checkpoint_path")
    sys.ρ = backend == :cuda ? ITensors.adapt(Array, rho) : rho
    if backend == :cuda
        sys.H0 = ITensors.adapt(Array, sys.H0)
        sys.VH = ITensors.adapt(Array, sys.VH)
        sys.VF = ITensors.adapt(Array, sys.VF)
    end
    write_mpo_checkpoint(sys, checkpoint_path)
end

open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "campaign" => string(campaign.name), "label" => string(spec.label),
        "task_index" => task_index, "diagnostic" => "fixed_initial_hamiltonian_sp2",
        "matrix_dimension" => N, "target_particles" => Ne,
        "spectral_lower" => bounds[1], "spectral_upper" => bounds[2],
        "campaign_spectral_lower" => campaign_bounds[1],
        "campaign_spectral_upper" => campaign_bounds[2],
        "backend" => string(backend), "device_name" => device_name,
        "device_total_memory_bytes" => device_total_memory,
        "device_free_memory_before_bytes" => device_free_memory_before,
        "device_free_memory_after_bytes" => device_free_memory_after,
        "itensors_tol" => params.itensors_tol, "itensors_maxdim" => maxdim,
        "purification_steps" => params.purification_steps,
        "trace_tolerance" => trace_tolerance,
        "idempotency_tolerance" => idempotency_tolerance,
        "hermiticity_tolerance" => hermiticity_tolerance,
        "initial_rho_max_chi" => maxlinkdim(initialization.value),
        "initial_rho_mean_chi" => mean_bond_dimension(initialization.value),
        "initialization_time_s" => initialization.time,
        "initialization_allocations_bytes" => initialization.bytes,
        "initialization_gc_time_s" => initialization.gctime,
        "converged" => converged, "termination_reason" => string(termination_reason),
        "iterations" => last_iteration, "final_rho_max_chi" => maxlinkdim(rho),
        "final_rho_mean_chi" => mean_bond_dimension(rho), "final_trace" => last_trace,
        "final_trace_error" => last_trace_error,
        "final_idempotency_residual" => last_idempotency,
        "final_hermiticity_residual" => last_hermiticity,
        "checkpoint_path" => isempty(checkpoint_path) ? "" : abspath(checkpoint_path),
    ))
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "started_at" => string(started_at), "finished_at" => string(now(UTC)),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(), "git_revision" => git_revision(repo_root),
        "project_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Project.toml")))),
        "manifest_sha1" => bytes2hex(sha1(read(joinpath(repo_root, "Manifest.toml")))),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "local"),
        "slurm_array_task_id" => get(ENV, "SLURM_ARRAY_TASK_ID", "local"),
    ))
end

println("Fixed-Hamiltonian SP2 diagnostic: label=$(spec.label) maxdim=$maxdim")
println("converged=$converged termination=$termination_reason iterations=$last_iteration")
println("final trace error=$last_trace_error idempotency=$last_idempotency chi=$(maxlinkdim(rho))")
println("Result directory: $output")
