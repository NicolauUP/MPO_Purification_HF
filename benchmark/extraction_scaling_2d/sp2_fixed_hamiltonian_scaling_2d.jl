#!/usr/bin/env julia

"""Fixed-Hamiltonian, large-square SP2 scalability and invariant diagnostic.

No SCF iteration and no dense matrix is constructed. The checkerboard
potential gives a deliberately gapped, noninteracting reference Hamiltonian.
For every requested square size and density-MPO cap this program records the
SP2 trace, idempotency, Hermiticity, branch, bond dimension, time, and
allocation after every polynomial. It then verifies the final fused-adjacency
Hartree field at a small set of directly measured sites.
"""

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using SHA

_csv(value) = "\"" * replace(string(value), '\"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))

function _bond_dimensions(mpo::MPO)
    length(mpo) <= 1 && return Int[]
    [dim(commonind(mpo[index], mpo[index + 1])) for index in 1:(length(mpo) - 1)]
end

_mean_chi(mpo::MPO) = begin
    dimensions = _bond_dimensions(mpo)
    isempty(dimensions) ? 1.0 : sum(dimensions) / length(dimensions)
end

function _parse_ints(value)
    Tuple(parse(Int, strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end

function parse_arguments(arguments)
    configuration = (
        output=nothing, side_levels=(6, 8, 10), maxdims=(256,),
        itensors_tol=1e-14, steps=50, padding=0.5,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-levels" && (configuration = merge(configuration, (side_levels=_parse_ints(arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdims" && (configuration = merge(configuration, (maxdims=_parse_ints(arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (configuration = merge(configuration, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (configuration = merge(configuration, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--padding" && (configuration = merge(configuration, (padding=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(configuration.output) && throw(ArgumentError("--output DIRECTORY is required"))
    all(level -> 1 <= level <= 10, configuration.side_levels) || throw(ArgumentError("side levels must lie in 1:10"))
    all(>(0), configuration.maxdims) || throw(ArgumentError("all maxdims must be positive"))
    isfinite(configuration.itensors_tol) && configuration.itensors_tol > 0 || throw(ArgumentError("itensors tolerance must be positive and finite"))
    configuration.steps > 0 || throw(ArgumentError("steps must be positive"))
    isfinite(configuration.padding) && configuration.padding >= 0 || throw(ArgumentError("padding must be finite and nonnegative"))
    configuration
end

function _parameters(side_level, maxdim, tolerance, steps)
    ParametersSquare(
        L=2side_level, t=(-0.6, -0.35), U=0.3,
        W=(x, y) -> iseven(Int(x) + Int(y)) ? 0.6 : -0.6,
        S=nothing, tci_tol=1e-10, itensors_tol=tolerance,
        itensors_maxdim=maxdim, density=0.5,
        purification_steps=steps, scf_mixing=0.5, scf_tol=0.1,
        scf_max_iterations=1,
    )
end

function _probe_coordinates(side)
    unique(((0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1),
        (div(side, 2), div(side, 2))))
end

function _direct_hartree(sys::System, site::Int)
    params = sys.params
    params.U * sum(
        real(MatrixChecker(sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states))
        for neighbour in values(square_neighbours(site, params.L)) if !isnothing(neighbour)
    )
end

function _hartree_probes(sys::System)
    field = extract_hartree_mpo_binary_carry_square_adjacency(sys)
    side = 2 ^ div(sys.params.L, 2)
    values, errors = String[], Float64[]
    for (x, y) in _probe_coordinates(side)
        site = square_lattice_index(x, y, sys.params.L)
        direct = _direct_hartree(sys, site)
        observed = real(MatrixChecker(field, sys.sites, site, site, sys.bra_states, sys.ket_states))
        push!(errors, abs(observed - direct))
        push!(values, @sprintf("(%d,%d):direct=%.16e;adjacency=%.16e;error=%+.3e",
            x, y, direct, observed, observed - direct))
    end
    return (
        max_error=maximum(errors), mean_error=sum(errors) / length(errors),
        values=join(values, ';'), field_max_chi=maxlinkdim(field), field_mean_chi=_mean_chi(field),
    )
end

function _write_metadata(output, configuration)
    project = Base.active_project()
    open(joinpath(output, "metadata.txt"), "w") do io
        println(io, "started_at_utc=", Dates.now(Dates.UTC))
        println(io, "active_project=", project)
        println(io, "project_sha256=", bytes2hex(sha256(read(project))))
        println(io, "julia_version=", VERSION)
        println(io, "cpu_name=", Sys.CPU_NAME)
        println(io, "julia_threads=", Threads.nthreads())
        println(io, "blas_threads=", BLAS.get_num_threads())
        println(io, "side_levels=", join(configuration.side_levels, ','))
        println(io, "maxdims=", join(configuration.maxdims, ','))
        println(io, "itensors_tol=", configuration.itensors_tol)
        println(io, "purification_steps=", configuration.steps)
        println(io, "potential=checkerboard_plus_minus_0.6")
        println(io, "hamiltonian=noninteracting_square_t=(-0.6,-0.35)")
        println(io, "base_spectral_bounds=(-3.5,3.5)")
        println(io, "spectral_padding=", configuration.padding)
        println(io, "spectral_bound_reason=2*(abs(tx)+abs(ty))+abs(W)=2.5 < 3.5 for open boundaries; padding is applied")
        println(io, "note=no_dense_matrix_or_scf; Hartree probes compare field to the final MPO density")
    end
end

function _run_case(iterations, summary, side_level, maxdim, configuration)
    params = _parameters(side_level, maxdim, configuration.itensors_tol, configuration.steps)
    sys = System(params)
    N = 2 ^ params.L
    Ne = round(Int, N * params.density)
    bounds = (-3.5 - configuration.padding, 3.5 + configuration.padding)
    rho = construct_rho_0(sys, params, bounds...; method=:sp2)
    trace_tolerance = MPO_MeanField._sp2_trace_tolerance(params, Ne)
    hermiticity_tolerance = MPO_MeanField._sp2_hermiticity_tolerance(params)
    idempotency_tolerance = 1e-3
    previous_idempotency, stagnant_steps = Inf, 0
    termination = :max_iterations
    converged = false

    for iteration in 1:params.purification_steps
        timing = @timed begin
            rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            trace_value = real(tr(rho))
            trace_squared = real(tr(rho_squared))
            idempotency = MPO_MeanField.idempotency_residual(rho, rho_squared)
            hermiticity = MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
            branch = MPO_MeanField._sp2_branch(trace_value, trace_squared, Ne, trace_tolerance)
            (rho_squared, trace_value, idempotency, hermiticity, branch)
        end
        rho_squared, trace_value, idempotency, hermiticity, branch = timing.value
        _row(iterations, (
            side_level, params.L, N, maxdim, iteration, branch,
            trace_value, abs(trace_value - Ne), idempotency, hermiticity,
            maxlinkdim(rho), maxlinkdim(rho_squared), _mean_chi(rho), _mean_chi(rho_squared),
            timing.time, timing.bytes, timing.gctime,
        ))
        flush(iterations)

        if abs(trace_value - Ne) <= trace_tolerance &&
           idempotency < idempotency_tolerance && hermiticity <= hermiticity_tolerance
            termination, converged = :idempotency_threshold, true
            break
        end
        if idempotency >= previous_idempotency * (1 - 1e-12)
            stagnant_steps += 1
        else
            stagnant_steps = 0
        end
        if stagnant_steps >= 8
            termination = :stagnation
            break
        end
        previous_idempotency = idempotency
        rho = branch == :square ? rho_squared : +(
            2.0 * rho, -rho_squared; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
    end
    sys.ρ = rho
    probes_timing = @timed _hartree_probes(sys)
    probes = probes_timing.value
    final_square = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    final_idempotency = MPO_MeanField.idempotency_residual(rho, final_square)
    final_hermiticity = MPO_MeanField._relative_mpo_residual(rho, ITensors.dag(rho), params)
    _row(summary, (
        side_level, params.L, N, Ne, maxdim, params.itensors_tol, bounds[1], bounds[2],
        converged, termination, maxlinkdim(rho), _mean_chi(rho),
        real(tr(rho)), abs(real(tr(rho)) - Ne), final_idempotency, final_hermiticity,
        probes_timing.time, probes_timing.bytes, probes.field_max_chi, probes.field_mean_chi,
        probes.max_error, probes.mean_error, probes.values,
    ))
    flush(summary)
end

function run_scaling(; output, side_levels, maxdims, itensors_tol, steps, padding)
    output = abspath(output)
    mkpath(output)
    configuration = (; output, side_levels, maxdims, itensors_tol, steps, padding)
    _write_metadata(output, configuration)
    open(joinpath(output, "iterations.csv"), "w") do iterations
        _row(iterations, (
            "L_side", "L_total", "N", "rho_maxdim", "iteration", "branch",
            "trace", "trace_error", "idempotency_residual", "hermiticity_residual",
            "rho_chi", "rho_squared_chi", "rho_mean_chi", "rho_squared_mean_chi",
            "step_time_s", "step_allocations_bytes", "step_gc_time_s",
        ))
        open(joinpath(output, "summary.csv"), "w") do summary
            _row(summary, (
                "L_side", "L_total", "N", "Ne", "rho_maxdim", "itensors_tol", "spectral_lower", "spectral_upper",
                "converged", "termination_reason", "final_rho_chi", "final_rho_mean_chi",
                "trace", "trace_error", "idempotency_residual", "hermiticity_residual",
                "hartree_time_s", "hartree_allocations_bytes", "hartree_field_max_chi", "hartree_field_mean_chi",
                "direct_probe_max_abs_error", "direct_probe_mean_abs_error", "direct_probe_values",
            ))
            open(joinpath(output, "errors.csv"), "w") do errors
                _row(errors, ("L_side", "L_total", "rho_maxdim", "stage", "error"))
                for side_level in side_levels, maxdim in maxdims
                    println("L_side=$side_level L_total=$(2side_level) maxdim=$maxdim: fixed-Hamiltonian SP2")
                    try
                        _run_case(iterations, summary, side_level, maxdim, configuration)
                    catch error
                        _row(errors, (side_level, 2side_level, maxdim, "case", sprint(showerror, error, catch_backtrace())))
                        flush(errors)
                        println(stderr, "L_side=$side_level maxdim=$maxdim failed; details recorded in errors.csv")
                    end
                end
            end
        end
    end
    println("Fixed-Hamiltonian SP2 scaling artefacts written to $output")
end

function print_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/sp2_fixed_hamiltonian_scaling_2d.jl --output DIRECTORY [options]")
    println("  --side-levels 6,8,10 --maxdims 256 --itensors-tol 1e-14 --steps 50 --padding 0.5")
    println("  no dense matrix and no SCF are constructed; this is a fixed-Hamiltonian resource/invariant ladder")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_arguments(ARGS)
    isnothing(configuration) ? print_usage() : run_scaling(; configuration...)
end
