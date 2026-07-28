#!/usr/bin/env julia

"""Bounded single-case 2D square Hartree--Fock SCF pilot.

This is intentionally a *single-run* diagnostic, not a parameter sweep. It
writes sufficient numerical history and provenance to assess whether a chosen
MPO cap, truncation tolerance, and initial seed support a trustworthy later
campaign. To avoid accidental runs beyond the current scaling study,
`side_level` is capped at ten (`side = 1024`, `N = 1_048_576`).
"""

using Dates
using ITensors, ITensorMPS
using MPO_MeanField
using Printf
using SHA

function _requested_cuda_backend(arguments)
    index = findfirst(==("--backend"), arguments)
    return !isnothing(index) && index < length(arguments) &&
        lowercase(arguments[index + 1]) == "cuda"
end

# Load CUDA before `run_pilot` is defined. Loading it with `@eval` from inside
# that compiled function creates a newer world-age binding on Julia 1.12.
const CUDA_BACKEND = if _requested_cuda_backend(ARGS)
    @eval import CUDA
    CUDA
else
    nothing
end

_csv(value) = "\"" * replace(string(value), '\"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))
_toml_string(value) = "\"" * replace(string(value), '\"' => "\\\"") * "\""
_optional(value) = isnothing(value) ? "" : string(value)

function parse_arguments(arguments)
    config = (
        output=nothing, side_level=6, tx=-0.6, ty=-0.35, U=0.3,
        density=0.5, potential=:checkerboard, potential_amplitude=0.6,
        seed=:uniform, seed_amplitude=0.05, tci_tol=1e-10,
        itensors_tol=1e-14, maxdim=256, steps=50, padding=0.5,
        scf_mixing=0.5, scf_tol=0.1, scf_max_iterations=30,
        backend=:cpu, sp2_idempotency_tolerance=1e-3,
        sp2_relative_trace_tolerance=1e-6, observables=:compact,
        square_fock_method=:binary_carry,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (config = merge(config, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (config = merge(config, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--tx" && (config = merge(config, (tx=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--ty" && (config = merge(config, (ty=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--U" && (config = merge(config, (U=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--density" && (config = merge(config, (density=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--potential" && (config = merge(config, (potential=Symbol(arguments[index + 1]),)); index += 2; continue)
        argument == "--potential-amplitude" && (config = merge(config, (potential_amplitude=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--seed" && (config = merge(config, (seed=Symbol(arguments[index + 1]),)); index += 2; continue)
        argument == "--seed-amplitude" && (config = merge(config, (seed_amplitude=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--tci-tol" && (config = merge(config, (tci_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (config = merge(config, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdim" && (config = merge(config, (maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (config = merge(config, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--padding" && (config = merge(config, (padding=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-mixing" && (config = merge(config, (scf_mixing=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-tol" && (config = merge(config, (scf_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--scf-max-iterations" && (config = merge(config, (scf_max_iterations=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--backend" && (config = merge(config, (backend=Symbol(lowercase(arguments[index + 1])),)); index += 2; continue)
        argument == "--sp2-idempotency-tolerance" && (config = merge(config, (sp2_idempotency_tolerance=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--sp2-relative-trace-tolerance" && (config = merge(config, (sp2_relative_trace_tolerance=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--observables" && (config = merge(config, (observables=Symbol(lowercase(arguments[index + 1])),)); index += 2; continue)
        argument == "--square-fock-method" && (config = merge(config, (square_fock_method=Symbol(lowercase(arguments[index + 1])),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(config.output) && throw(ArgumentError("--output DIRECTORY is required"))
    1 <= config.side_level <= 10 || throw(ArgumentError("side level must lie in 1:10 for this bounded SCF pilot"))
    config.potential in (:none, :checkerboard) || throw(ArgumentError("--potential must be none or checkerboard"))
    config.seed in (:uniform, :checkerboard_plus, :checkerboard_minus) || throw(ArgumentError("--seed must be uniform, checkerboard_plus, or checkerboard_minus"))
    config.backend in (:cpu, :cuda) || throw(ArgumentError("--backend must be cpu or cuda"))
    config.observables in (:compact, :full) ||
        throw(ArgumentError("--observables must be compact or full"))
    config.square_fock_method in (:tci, :binary_carry) || throw(ArgumentError(
        "--square-fock-method must be tci or binary_carry",
    ))
    all(isfinite, (config.tx, config.ty, config.U, config.density, config.potential_amplitude, config.seed_amplitude, config.tci_tol, config.itensors_tol, config.padding, config.scf_mixing, config.scf_tol)) || throw(ArgumentError("all floating-point inputs must be finite"))
    config.potential_amplitude >= 0 || throw(ArgumentError("potential amplitude must be nonnegative"))
    config.seed_amplitude >= 0 || throw(ArgumentError("seed amplitude must be nonnegative"))
    config.tci_tol > 0 && config.itensors_tol > 0 || throw(ArgumentError("TCI and ITensor tolerances must be positive"))
    config.maxdim > 0 && config.steps > 0 && config.scf_max_iterations > 0 || throw(ArgumentError("maxdim, steps, and SCF iterations must be positive"))
    0 < config.density < 1 || throw(ArgumentError("density must lie strictly between zero and one"))
    0 <= config.scf_mixing <= 1 && config.scf_tol > 0 && config.padding >= 0 || throw(ArgumentError("invalid SCF mixing, tolerance, or padding"))
    config.sp2_idempotency_tolerance > 0 ||
        throw(ArgumentError("SP2 idempotency tolerance must be positive"))
    config.sp2_relative_trace_tolerance > 0 ||
        throw(ArgumentError("SP2 relative trace tolerance must be positive"))
    return config
end

function _checkerboard(amplitude, sign=1.0)
    (x, y) -> sign * (iseven(Int(x) + Int(y)) ? amplitude : -amplitude)
end

function _parameters(config)
    potential = config.potential == :none ? nothing : _checkerboard(config.potential_amplitude)
    seed = if config.seed == :uniform
        nothing
    elseif config.seed == :checkerboard_plus
        _checkerboard(config.seed_amplitude)
    else
        _checkerboard(config.seed_amplitude, -1.0)
    end
    return ParametersSquare(
        L=2config.side_level, t=(config.tx, config.ty), U=config.U,
        W=potential, S=seed, tci_tol=config.tci_tol,
        itensors_tol=config.itensors_tol, itensors_maxdim=config.maxdim,
        density=config.density, purification_steps=config.steps,
        scf_mixing=config.scf_mixing, scf_tol=config.scf_tol,
        scf_max_iterations=config.scf_max_iterations,
    )
end

function _write_metadata(path, config, params, bounds)
    open(path, "w") do io
        println(io, "format_version = 1")
        println(io, "created_at_utc = ", _toml_string(Dates.format(now(UTC), DateFormat("yyyy-mm-ddTHH:MM:SS"))))
        println(io, "runner = ", _toml_string("benchmark/scf_pilot_2d/scf_pilot_2d.jl"))
        println(io, "project_toml_sha256 = ", _toml_string(bytes2hex(sha256(read("Project.toml")))))
        println(io, "side_level = $(config.side_level)")
        println(io, "L_total = $(params.L)")
        println(io, "N = $(2^params.L)")
        println(io, "target_particles = $(round(Int, config.density * 2^params.L))")
        for name in (:tx, :ty, :U, :density, :potential_amplitude, :seed_amplitude, :tci_tol, :itensors_tol, :padding, :scf_mixing, :scf_tol)
            println(io, "$(name) = $(getfield(config, name))")
        end
        println(io, "potential = ", _toml_string(config.potential))
        println(io, "seed = ", _toml_string(config.seed))
        println(io, "itensors_maxdim = $(config.maxdim)")
        println(io, "purification_steps = $(config.steps)")
        println(io, "scf_max_iterations = $(config.scf_max_iterations)")
        println(io, "backend = ", _toml_string(config.backend))
        println(io, "sp2_idempotency_tolerance = $(config.sp2_idempotency_tolerance)")
        println(io, "sp2_relative_trace_tolerance = $(config.sp2_relative_trace_tolerance)")
        println(io, "observables = ", _toml_string(config.observables))
        println(io, "square_fock_method = ", _toml_string(config.square_fock_method))
        println(io, "spectral_lower = $(bounds[1])")
        println(io, "spectral_upper = $(bounds[2])")
    end
end

function _write_history(path, diagnostics)
    open(path, "w") do io
        _row(io, ("iteration", "trace", "vh_residual", "vf_residual", "rho_residual", "commutator_residual", "two_cycle_residual", "purification_converged", "purification_termination_reason", "purification_iterations", "rho_bond_dimension", "hartree_bond_dimension", "fock_bond_dimension", "effective_hamiltonian_bond_dimension", "energy_total"))
        for record in diagnostics.history
            _row(io, (record.iteration, record.trace, record.vh_residual, record.vf_residual, record.rho_residual, record.commutator_residual, record.two_cycle_residual, record.purification_converged, record.purification_termination_reason, record.purification_iterations, _optional(record.rho_bond_dimension), _optional(record.hartree_bond_dimension), _optional(record.fock_bond_dimension), _optional(record.effective_hamiltonian_bond_dimension), _optional(record.energy_total)))
        end
    end
end

function _write_post_scf_timings(path, timings)
    open(path, "w") do io
        for name in propertynames(timings)
            println(io, "$(name) = $(getfield(timings, name))")
        end
    end
end

function _probe_sites(side, L)
    coordinates = unique(((0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1), (div(side, 2), div(side, 2))))
    [(x, y, MPO_MeanField.square_lattice_index(x, y, L)) for (x, y) in coordinates]
end

function _compact_density_summary(sys)
    params = sys.params
    N = 2^params.L
    checkerboard = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = MPO_MeanField.square_lattice_decoder(Int(z), params.L)
            iseven(x + y) ? 1.0 : -1.0
        end,
        Float64,
        sys.sites,
        params.tci_tol,
    )
    weighted_density = apply(
        checkerboard, sys.ρ;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return (
        particle_number=real(tr(sys.ρ)),
        checkerboard_order=real(tr(weighted_density)) / N,
        density_min=nothing,
        density_max=nothing,
    )
end

function _full_density_summary(sys)
    params = sys.params
    N = 2^params.L
    density_sum = 0.0
    density_min, density_max = Inf, -Inf
    checkerboard_sum = 0.0
    for site in 1:N
        value = real(MatrixChecker(sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states))
        x, y = MPO_MeanField.square_lattice_decoder(site - 1, params.L)
        density_sum += value
        density_min = min(density_min, value)
        density_max = max(density_max, value)
        checkerboard_sum += iseven(x + y) ? value : -value
    end
    return (
        particle_number=density_sum,
        checkerboard_order=checkerboard_sum / N,
        density_min=density_min,
        density_max=density_max,
    )
end

function _compact_energy(sys)
    params = sys.params
    final_hartree, final_fock, _ = MPO_MeanField._extract_mean_fields_with_components(sys)
    H0rho = apply(
        sys.H0, sys.ρ;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    VHrho = apply(
        final_hartree, sys.ρ;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    VFrho = apply(
        final_fock, sys.ρ;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    kinetic = real(tr(H0rho))
    hartree = real(tr(VHrho)) / 2
    fock = real(tr(VFrho)) / 2
    interaction = hartree + fock
    return (
        kinetic=kinetic,
        hartree=hartree,
        fock=fock,
        interaction=interaction,
        total=kinetic + interaction,
    )
end

function _measure_probes(sys, side, L)
    [(
        x=x,
        y=y,
        site=site,
        density=real(MatrixChecker(
            sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states,
        )),
        stored_hartree=real(MatrixChecker(
            sys.VH, sys.sites, site, site, sys.bra_states, sys.ket_states,
        )),
    ) for (x, y, site) in _probe_sites(side, L)]
end

function _write_observables(path, sys, converged, bounds; mode=:compact)
    params = sys.params
    N = 2^params.L
    side = 2^div(params.L, 2)
    density_timing = @timed(
        mode == :compact ? _compact_density_summary(sys) : _full_density_summary(sys)
    )
    density = density_timing.value
    energy_timing = @timed(
        mode == :compact ? _compact_energy(sys) : nearest_neighbor_hf_energy_square(sys)
    )
    energy = energy_timing.value
    effective_timing = @timed(
        +(sys.H0, sys.VH, sys.VF;
            cutoff=params.itensors_tol,
            maxdim=params.itensors_maxdim,
        )
    )
    effective_hamiltonian = effective_timing.value
    idempotency_timing = @timed MPO_MeanField._mpo_relative_change(
        apply(sys.ρ, sys.ρ; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim),
        sys.ρ,
        params,
    )
    idempotency = idempotency_timing.value
    stationarity_timing = @timed MPO_MeanField._scf_commutator_residual(
        effective_hamiltonian, sys.ρ, params,
    )
    stationarity = stationarity_timing.value
    probes_timing = @timed _measure_probes(sys, side, params.L)
    probes = probes_timing.value
    open(path, "w") do io
        println(io, "observable_mode = ", _toml_string(mode))
        println(io, "energy_evaluation = ", _toml_string(
            mode == :compact ? "fresh_mean_field_trace" : "explicit_nearest_neighbor",
        ))
        println(io, "scf_converged = $(converged)")
        println(io, "scf_termination_reason = ", _toml_string(scf_diagnostics(sys).termination_reason))
        println(io, "spectral_lower = $(bounds[1])")
        println(io, "spectral_upper = $(bounds[2])")
        println(io, "particle_number = $(density.particle_number)")
        println(io, "particle_number_error = $(density.particle_number - round(Int, params.density * N))")
        println(io, "density_extrema_computed = $(mode == :full)")
        if mode == :full
            println(io, "density_min = $(density.density_min)")
            println(io, "density_max = $(density.density_max)")
        end
        println(io, "checkerboard_density_order = $(density.checkerboard_order)")
        println(io, "idempotency_residual = $(idempotency)")
        println(io, "stationarity_residual = $(stationarity)")
        println(io, "rho_bond_dimension = $(maxlinkdim(sys.ρ))")
        println(io, "hartree_bond_dimension = $(maxlinkdim(sys.VH))")
        println(io, "fock_bond_dimension = $(maxlinkdim(sys.VF))")
        println(io, "effective_hamiltonian_bond_dimension = $(maxlinkdim(effective_hamiltonian))")
        for name in (:kinetic, :hartree, :fock, :interaction, :total)
            println(io, "energy_$(name) = $(getfield(energy, name))")
        end
        for probe in probes
            println(io, "\n[[density_probes]]")
            println(io, "x = $(probe.x)\ny = $(probe.y)\nsite = $(probe.site)")
            println(io, "density = $(probe.density)")
            println(io, "stored_hartree = $(probe.stored_hartree)")
        end
    end
    return (
        density_summary_time_s=density_timing.time,
        energy_time_s=energy_timing.time,
        effective_hamiltonian_time_s=effective_timing.time,
        idempotency_time_s=idempotency_timing.time,
        stationarity_time_s=stationarity_timing.time,
        probes_time_s=probes_timing.time,
    )
end

function run_pilot(config)
    output = abspath(config.output)
    mkpath(output)
    params = _parameters(config)
    potential_bounds = config.potential == :none ? (0.0, 0.0) : (-config.potential_amplitude, config.potential_amplitude)
    bounds = square_scf_spectral_bounds(params; potential_bounds=potential_bounds, margin=config.padding)
    _write_metadata(joinpath(output, "metadata.toml"), config, params, bounds)
    if config.backend == :cuda
        isnothing(CUDA_BACKEND) &&
            error("CUDA backend must be selected with the command-line option --backend cuda")
        CUDA_BACKEND.functional() || error("CUDA is not functional on this node")
    end
    to_device = config.backend == :cuda ? CUDA_BACKEND.cu : identity
    to_host = config.backend == :cuda ? ITensors.cpu : identity
    synchronize_backend = config.backend == :cuda ? CUDA_BACKEND.synchronize : (() -> nothing)
    gpu_total_memory = config.backend == :cuda ? CUDA_BACKEND.total_memory() : 0
    sys = System(params)
    phase_path = joinpath(output, "phase_timings.csv")
    converged = open(joinpath(output, "progress.txt"), "w") do progress
        open(phase_path, "w") do phases
            _row(phases, (
                "iteration", "initialization_time_s", "purification_time_s",
                "density_to_host_time_s", "mean_field_time_s",
                "fields_to_device_time_s", "device_diagnostics_time_s",
                "residuals_time_s", "mixing_time_s",
                "measured_iteration_time_s", "purification_iterations",
                "rho_bond_dimension", "hartree_bond_dimension",
                "fock_bond_dimension", "effective_hamiltonian_bond_dimension",
                "gpu_free_memory_bytes", "gpu_used_memory_bytes",
                "gpu_total_memory_bytes",
            ))
            phase_callback = record -> begin
                free_memory = config.backend == :cuda ? CUDA_BACKEND.free_memory() : 0
                _row(phases, (
                    record.iteration, record.initialization_time_s,
                    record.purification_time_s, record.density_to_host_time_s,
                    record.mean_field_time_s, record.fields_to_device_time_s,
                    record.device_diagnostics_time_s,
                    record.residuals_time_s, record.mixing_time_s,
                    record.measured_iteration_time_s,
                    record.purification_iterations, record.rho_bond_dimension,
                    record.hartree_bond_dimension, record.fock_bond_dimension,
                    record.effective_hamiltonian_bond_dimension, free_memory,
                    gpu_total_memory - free_memory, gpu_total_memory,
                ))
                flush(phases)
            end
            run_scf!(
                sys, bounds...;
                purification_method=:sp2,
                square_fock_method=config.square_fock_method,
                sp2_idempotency_tolerance=config.sp2_idempotency_tolerance,
                sp2_trace_tolerance=config.sp2_relative_trace_tolerance *
                    round(Int, params.density * 2^params.L),
                verbose=:all,
                io=progress,
                overwrite_progress=false,
                record_energy=config.backend == :cpu,
                to_gpu=to_device,
                to_cpu=to_host,
                phase_callback=phase_callback,
                phase_synchronize=synchronize_backend,
            )
        end
    end
    post_scf_started = time_ns()
    cuda_status_time = if config.backend == :cuda
        @elapsed begin
            synchronize_backend()
            open(joinpath(output, "cuda_status.txt"), "w") do io
                CUDA_BACKEND.versioninfo(io)
                println(io)
                CUDA_BACKEND.pool_status(io)
            end
        end
    else
        0.0
    end
    fields_to_host_time = if config.backend == :cuda
        @elapsed begin
            sys.H0 = to_host(sys.H0)
            sys.VH = to_host(sys.VH)
            sys.VF = to_host(sys.VF)
            synchronize_backend()
        end
    else
        0.0
    end
    history_time = @elapsed(
        _write_history(joinpath(output, "scf_history.csv"), scf_diagnostics(sys))
    )
    observable_timings = _write_observables(
        joinpath(output, "observables.toml"), sys, converged, bounds;
        mode=config.observables,
    )
    post_scf_total_time = (time_ns() - post_scf_started) / 1e9
    _write_post_scf_timings(
        joinpath(output, "post_scf_timings.toml"),
        merge((
            cuda_status_time_s=cuda_status_time,
            fields_to_host_time_s=fields_to_host_time,
            history_time_s=history_time,
        ), observable_timings, (total_time_s=post_scf_total_time,)),
    )
    println("SCF pilot completed: converged=$converged termination=$(scf_diagnostics(sys).termination_reason)")
    println("Output directory: $output")
    return converged
end

function main(arguments)
    config = parse_arguments(arguments)
    if isnothing(config)
        println("Usage: scf_pilot_2d.jl --output DIRECTORY [options]")
        println("Options: --side-level 1..10 --potential none|checkerboard --seed uniform|checkerboard_plus|checkerboard_minus --backend cpu|cuda --observables compact|full --square-fock-method tci|binary_carry")
        return
    end
    run_pilot(config)
end

main(ARGS)
