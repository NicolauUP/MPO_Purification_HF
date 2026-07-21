#!/usr/bin/env julia

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using SHA

const DEFAULT_SIDE_LEVELS = Tuple(2:10)
const DEFAULT_SOURCES = (:smooth, :checkerboard_exact)

_csv(value) = "\"" * replace(string(value), '"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))
_median(values) = sort(values)[cld(length(values), 2)]

function _bond_dimensions(mpo::MPO)
    length(mpo) <= 1 && return Int[]
    [dim(commonind(mpo[i], mpo[i + 1])) for i in 1:(length(mpo) - 1)]
end
_mean_chi(mpo::MPO) = isempty(_bond_dimensions(mpo)) ? 1.0 : sum(_bond_dimensions(mpo)) / length(_bond_dimensions(mpo))

function _relative_difference(left::MPO, right::MPO, params)
    difference = +(left, -right; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    numerator = sqrt(max(0.0, real(inner(difference, difference))))
    denominator = max(sqrt(max(0.0, real(inner(right, right)))), sqrt(eps(Float64)))
    numerator / denominator
end

_hermiticity(field::MPO, params) = _relative_difference(field, ITensors.dag(field), params)

function _params(total_bits::Int; staggered::Bool=false, itensors_tol::Float64=1e-12)
    ParametersSquare(
        L=total_bits, t=(-0.6, -0.35), U=0.3,
        W=staggered ? ((x, y) -> iseven(Int(x) + Int(y)) ? 0.6 : -0.6) : nothing,
        S=nothing, tci_tol=1e-10, itensors_tol=itensors_tol, itensors_maxdim=128,
        density=0.5, purification_steps=50, scf_mixing=0.5, scf_tol=0.1,
        scf_max_iterations=5,
    )
end

function _smooth_system(total_bits::Int; itensors_tol::Float64=1e-12)
    params = _params(total_bits; itensors_tol=itensors_tol)
    side = 2 ^ div(total_bits, 2)
    sys = System(params)
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            0.45 + 0.08 * cospi(x / (side - 1)) + 0.05 * sinpi(y / (side - 1))
        end, Float64, sys.sites, params.tci_tol,
    )
    horizontal = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            x < side - 1 ? 0.03 * sinpi((x + 2y) / (2side - 2)) : 0.0
        end, Float64, sys.sites, params.tci_tol,
    )
    vertical = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            y < side - 1 ? 0.025 * cospi((2x + y) / (2side - 2)) : 0.0
        end, Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L, T_U, T_D = sys.translations
    rho_right = apply(horizontal, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_left = apply(T_L, ITensors.dag(horizontal); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_up = apply(vertical, T_U; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_down = apply(T_D, ITensors.dag(vertical); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.ρ = +(diagonal, rho_right, rho_left, rho_up, rho_down;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return sys, "deterministic_smooth_synthetic", "not_applicable"
end

"""Exact, bounded checkerboard CDW density used to isolate Hartree extraction.

This is deliberately a density input, not the result of an SP2 or SCF solve:
the benchmark measures the field kernels even when current square-MPO SP2
stagnates. It has the physical bounds `0 <= n(x,y) <= 1` and a compact QTT
description at every benchmark size.
"""
function _checkerboard_exact_system(total_bits::Int; itensors_tol::Float64=1e-12)
    params = _params(total_bits; itensors_tol=itensors_tol)
    delta = 0.2
    sys = System(params)
    sys.ρ = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), total_bits)
            0.5 + (iseven(x + y) ? delta : -delta)
        end,
        Float64,
        sys.sites,
        params.tci_tol,
    )
    return sys, "exact_checkerboard_cdw", @sprintf("n=0.5%+.3f*(-1)^(x+y)", delta)
end

function _sp2_system(total_bits::Int; itensors_tol::Float64=1e-12)
    params = _params(total_bits; staggered=true, itensors_tol=itensors_tol)
    sys = System(params)
    bounds = (-3.5, 3.5)
    initial = construct_rho_0(sys, params, bounds...; method=:sp2)
    result = perform_purification(initial, params; method=:sp2, verbose=0, spectral_bounds=bounds)
    result.converged || error("SP2 preparation did not converge: $(result.termination_reason), trace_error=$(result.trace_error)")
    sys.ρ = result.rho
    details = @sprintf("sp2:iterations=%d;trace_error=%.3e;idempotency=%.3e;max_chi=%d",
        result.iterations, something(result.trace_error, NaN), result.idempotency_residual,
        result.final_bond_dimension)
    return sys, "noninteracting_gapped_sp2", details
end

function _system(source::Symbol, total_bits::Int; itensors_tol::Float64=1e-12)
    source == :smooth && return _smooth_system(total_bits; itensors_tol=itensors_tol)
    source == :checkerboard_exact && return _checkerboard_exact_system(total_bits; itensors_tol=itensors_tol)
    source == :sp2_gapped && return _sp2_system(total_bits; itensors_tol=itensors_tol)
    throw(ArgumentError("unknown source $source"))
end

const KERNELS = (
    (:hartree_four_carry, extract_hartree_mpo_binary_carry_square),
    (:hartree_adjacency, extract_hartree_mpo_binary_carry_square_adjacency),
)

function _probe_coordinates(side::Int)
    unique(((0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1),
        (div(side, 2), div(side, 2))))
end

function _direct_hartree(sys::System, site::Int)
    params = sys.params
    sum(
        real(MatrixChecker(sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states))
        for neighbour in values(square_neighbours(site, params.L))
        if !isnothing(neighbour)
    ) * params.U
end

function _direct_probe_diagnostics(field::MPO, sys::System, kernel::Symbol)
    L = sys.params.L
    side = 2 ^ div(L, 2)
    coordinates = _probe_coordinates(side)
    values = String[]
    errors = Float64[]
    for (x, y) in coordinates
        site = square_lattice_index(x, y, L)
        direct = _direct_hartree(sys, site)
        observed = real(MatrixChecker(field, sys.sites, site, site, sys.bra_states, sys.ket_states))
        push!(errors, abs(observed - direct))
        push!(values, @sprintf("(%d,%d):direct=%.16e;%s=%.16e;error=%+.3e",
            x, y, direct, kernel, observed, observed - direct))
    end
    return maximum(errors), sum(errors) / length(errors), join(values, ';')
end

function _time_kernel(kernel, sys, warmups, repetitions)
    for _ in 1:warmups
        kernel(sys)
    end
    results = NamedTuple[]
    for repetition in 1:repetitions
        GC.gc(true)
        timing = @timed kernel(sys)
        field = timing.value
        push!(results, (repetition=repetition, time_s=timing.time, bytes=timing.bytes,
            gc_time_s=timing.gctime, max_chi=maxlinkdim(field), mean_chi=_mean_chi(field),
            hermiticity=_hermiticity(field, sys.params), field=field))
    end
    results
end

function _parse_symbols(value)
    Tuple(Symbol(strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end
function _parse_levels(value)
    Tuple(parse(Int, strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end

function parse_arguments(arguments)
    configuration = (output=nothing, side_levels=DEFAULT_SIDE_LEVELS, sources=DEFAULT_SOURCES,
        itensors_tol=1e-12, warmups=1, repetitions_small=5, repetitions_large=3)
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-levels" && (configuration = merge(configuration, (side_levels=_parse_levels(arguments[index + 1]),)); index += 2; continue)
        argument == "--sources" && (configuration = merge(configuration, (sources=_parse_symbols(arguments[index + 1]),)); index += 2; continue)
        argument == "--itensors-tol" && (configuration = merge(configuration, (itensors_tol=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--warmups" && (configuration = merge(configuration, (warmups=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--repetitions-small" && (configuration = merge(configuration, (repetitions_small=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--repetitions-large" && (configuration = merge(configuration, (repetitions_large=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(configuration.output) && throw(ArgumentError("--output DIRECTORY is required"))
    all(level -> 1 <= level <= 10, configuration.side_levels) || throw(ArgumentError("side levels must lie in 1:10"))
    isfinite(configuration.itensors_tol) && configuration.itensors_tol > 0 || throw(ArgumentError("itensors tolerance must be positive and finite"))
    configuration
end

function _metadata(output, configuration)
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
        println(io, "total_bits=", join((2level for level in configuration.side_levels), ','))
        println(io, "sources=", join(string.(configuration.sources), ','))
        println(io, "itensors_tol=", configuration.itensors_tol)
        println(io, "note=peak RSS is process-global; see process_time.txt")
    end
end

function run_benchmark(; output, side_levels, sources, itensors_tol, warmups, repetitions_small, repetitions_large)
    output = abspath(output)
    mkpath(output)
    configuration = (; output, side_levels, sources, itensors_tol, warmups, repetitions_small, repetitions_large)
    _metadata(output, configuration)
    open(joinpath(output, "samples.csv"), "w") do samples
        _row(samples, ("source", "L_side", "L_total", "N", "kernel", "repetition", "time_s", "allocations_bytes", "gc_time_s", "rho_max_chi", "rho_mean_chi", "field_max_chi", "field_mean_chi", "hermiticity_residual"))
        open(joinpath(output, "summary.csv"), "w") do summary
            _row(summary, ("source", "source_details", "L_side", "L_total", "N", "kernel", "repetitions", "median_time_s", "minimum_time_s", "median_allocations_bytes", "median_gc_time_s", "rho_max_chi", "rho_mean_chi", "median_field_max_chi", "median_field_mean_chi", "hermiticity_residual", "direct_probe_max_abs_error", "direct_probe_mean_abs_error", "direct_probe_values"))
            open(joinpath(output, "errors.csv"), "w") do errors
                _row(errors, ("source", "L_side", "L_total", "stage", "error"))
                for source in sources, side_level in side_levels
                    total_bits = 2side_level
                    N = 2^total_bits
                    repetitions = side_level <= 5 ? repetitions_small : repetitions_large
                    println("source=$source L_side=$side_level L_total=$total_bits N=$N: preparing density")
                    try
                        sys, label, details = _system(source, total_bits; itensors_tol=itensors_tol)
                        results = Dict{Symbol,Vector{NamedTuple}}()
                        for (name, kernel) in KERNELS
                            println("source=$source L_side=$side_level kernel=$name: timing $repetitions repetitions")
                            measurements = _time_kernel(kernel, sys, warmups, repetitions)
                            results[name] = measurements
                            for item in measurements
                                _row(samples, (label, side_level, total_bits, N, name, item.repetition, item.time_s, item.bytes, item.gc_time_s, maxlinkdim(sys.ρ), _mean_chi(sys.ρ), item.max_chi, item.mean_chi, item.hermiticity))
                            end
                        end
                        for (name, _) in KERNELS
                            measurements = results[name]
                            direct_max_error, direct_mean_error, direct_values = _direct_probe_diagnostics(measurements[end].field, sys, name)
                            _row(summary, (label, details, side_level, total_bits, N, name, repetitions,
                                _median([item.time_s for item in measurements]), minimum(item.time_s for item in measurements),
                                _median([item.bytes for item in measurements]), _median([item.gc_time_s for item in measurements]),
                                maxlinkdim(sys.ρ), _mean_chi(sys.ρ), _median([item.max_chi for item in measurements]),
                                _median([item.mean_chi for item in measurements]), _median([item.hermiticity for item in measurements]),
                                direct_max_error, direct_mean_error, direct_values))
                        end
                        flush(samples); flush(summary)
                    catch error
                        _row(errors, (source, side_level, total_bits, "case", sprint(showerror, error, catch_backtrace())))
                        flush(errors)
                        println(stderr, "source=$source L_side=$side_level failed; details recorded in errors.csv")
                    end
                end
            end
        end
    end
    println("Benchmark artefacts written to $output")
end

function print_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/extraction_scaling_2d.jl --output DIRECTORY [options]")
    println("  --side-levels 2,3,...,10  (side=2^L_side; total bits=2L_side)")
    println("  --sources smooth,checkerboard_exact,sp2_gapped")
    println("  --itensors-tol 1e-12  (MPO cutoff; use 1e-14 for the SP2 diagnostic)")
    println("  --warmups N --repetitions-small N --repetitions-large N")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_arguments(ARGS)
    isnothing(configuration) ? print_usage() : run_benchmark(; configuration...)
end
