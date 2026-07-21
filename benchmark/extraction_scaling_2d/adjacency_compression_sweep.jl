#!/usr/bin/env julia

"""Explicit accuracy/rank sweep for the fused square Hartree adjacency field.

The reference is built without truncating `U * A_nn * n`. Each subsequent row
compresses a copy of that reference with one explicit `(cutoff, maxdim)`
policy. The resulting CSV is intended for plotting the accuracy--memory--time
Pareto frontier; it is not an SCF calculation.
"""

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using SHA

const DEFAULT_SIDE_LEVELS = (4, 6, 8, 10)
const DEFAULT_CUTOFFS = (1e-14, 1e-12, 1e-10, 1e-8)
const DEFAULT_MAXDIMS = (16, 32, 64, 128, 256)

_csv(value) = "\"" * replace(string(value), '\"' => "\"\"") * "\""
_row(io, values) = println(io, join(_csv.(values), ','))
_median(values) = sort(values)[cld(length(values), 2)]

function _bond_dimensions(mpo::MPO)
    length(mpo) <= 1 && return Int[]
    [dim(commonind(mpo[i], mpo[i + 1])) for i in 1:(length(mpo) - 1)]
end
_mean_chi(mpo::MPO) = isempty(_bond_dimensions(mpo)) ? 1.0 : sum(_bond_dimensions(mpo)) / length(_bond_dimensions(mpo))

function _params(total_bits::Int)
    ParametersSquare(
        L=total_bits, t=(-0.6, -0.35), U=0.3, W=nothing, S=nothing,
        tci_tol=1e-10, itensors_tol=1e-12, itensors_maxdim=128,
        density=0.5, purification_steps=50, scf_mixing=0.5, scf_tol=0.1,
        scf_max_iterations=5,
    )
end

function _smooth_system(total_bits::Int)
    params = _params(total_bits)
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
    return sys
end

function _probe_coordinates(side::Int)
    unique(((0, 0), (side - 1, 0), (0, side - 1), (side - 1, side - 1),
        (div(side, 2), div(side, 2))))
end

function _direct_hartree(sys::System, site::Int)
    sum(
        real(MatrixChecker(sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states))
        for neighbour in values(square_neighbours(site, sys.params.L))
        if !isnothing(neighbour)
    ) * sys.params.U
end

function _probe_diagnostics(field::MPO, sys::System)
    side = 2 ^ div(sys.params.L, 2)
    errors = Float64[]
    imaginary_parts = Float64[]
    values = String[]
    for (x, y) in _probe_coordinates(side)
        site = square_lattice_index(x, y, sys.params.L)
        direct = _direct_hartree(sys, site)
        coefficient = MatrixChecker(field, sys.sites, site, site, sys.bra_states, sys.ket_states)
        observed = real(coefficient)
        push!(errors, abs(observed - direct))
        push!(imaginary_parts, abs(imag(coefficient)))
        push!(values, @sprintf("(%d,%d):direct=%.16e;field=%.16e;error=%+.3e",
            x, y, direct, observed, observed - direct))
    end
    return maximum(imaginary_parts), maximum(errors), sum(errors) / length(errors), join(values, ';')
end

function _parse_levels(value)
    Tuple(parse(Int, strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end
function _parse_cutoffs(value)
    Tuple(parse(Float64, strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end

function parse_arguments(arguments)
    configuration = (
        output=nothing, side_levels=DEFAULT_SIDE_LEVELS, cutoffs=DEFAULT_CUTOFFS,
        maxdims=DEFAULT_MAXDIMS, repetitions=3,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (configuration = merge(configuration, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-levels" && (configuration = merge(configuration, (side_levels=_parse_levels(arguments[index + 1]),)); index += 2; continue)
        argument == "--cutoffs" && (configuration = merge(configuration, (cutoffs=_parse_cutoffs(arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdims" && (configuration = merge(configuration, (maxdims=_parse_levels(arguments[index + 1]),)); index += 2; continue)
        argument == "--repetitions" && (configuration = merge(configuration, (repetitions=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(configuration.output) && throw(ArgumentError("--output DIRECTORY is required"))
    all(level -> 1 <= level <= 10, configuration.side_levels) || throw(ArgumentError("side levels must lie in 1:10"))
    all(>(0), configuration.maxdims) || throw(ArgumentError("maxdims must be positive"))
    all(>=(0), configuration.cutoffs) || throw(ArgumentError("cutoffs must be nonnegative"))
    configuration.repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    configuration
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
        println(io, "side_levels=", join(configuration.side_levels, ','))
        println(io, "cutoffs=", join(configuration.cutoffs, ','))
        println(io, "maxdims=", join(configuration.maxdims, ','))
        println(io, "repetitions=", configuration.repetitions)
        println(io, "source=deterministic_smooth_synthetic")
        println(io, "reference=untruncated fused U*A_nn*n field")
    end
end

function _timed_compression(reference::MPO, cutoff::Float64, maxdim::Int, repetitions::Int)
    # Compile and exercise the exact requested policy outside timing. Different
    # max dimensions can dispatch through distinct code paths.
    truncate(copy(reference); cutoff=cutoff, maxdim=maxdim)
    samples = NamedTuple[]
    for _ in 1:repetitions
        GC.gc(true)
        timing = @timed truncate(copy(reference); cutoff=cutoff, maxdim=maxdim)
        push!(samples, (field=timing.value, time_s=timing.time, bytes=timing.bytes, gc_time_s=timing.gctime))
    end
    samples
end

function run_sweep(; output, side_levels, cutoffs, maxdims, repetitions)
    output = abspath(output)
    mkpath(output)
    configuration = (; output, side_levels, cutoffs, maxdims, repetitions)
    _write_metadata(output, configuration)
    open(joinpath(output, "summary.csv"), "w") do summary
        _row(summary, (
            "L_side", "L_total", "N", "policy", "cutoff", "maxdim", "repetitions",
            "median_time_s", "minimum_time_s", "median_allocations_bytes", "median_gc_time_s",
            "rho_max_chi", "rho_mean_chi", "field_max_chi", "field_mean_chi",
            "direct_probe_max_abs_imag", "direct_probe_max_abs_error", "direct_probe_mean_abs_error",
            "direct_probe_values",
        ))
        open(joinpath(output, "errors.csv"), "w") do errors
            _row(errors, ("L_side", "L_total", "stage", "error"))
            for side_level in side_levels
                total_bits = 2side_level
                N = 2^total_bits
                println("L_side=$side_level L_total=$total_bits N=$N: preparing smooth density")
                try
                    sys = _smooth_system(total_bits)
                    println("L_side=$side_level: constructing untruncated adjacency reference")
                    # Density preparation and compilation are intentionally
                    # outside the measured field-construction time.
                    extract_hartree_mpo_binary_carry_square_adjacency(sys)
                    GC.gc(true)
                    reference_timing = @timed extract_hartree_mpo_binary_carry_square_adjacency(sys)
                    reference = reference_timing.value
                    reference_error = _probe_diagnostics(reference, sys)
                    _row(summary, (
                        side_level, total_bits, N, "reference_untruncated", 0.0, "unbounded", 1,
                        reference_timing.time, reference_timing.time, reference_timing.bytes, reference_timing.gctime,
                        maxlinkdim(sys.ρ), _mean_chi(sys.ρ), maxlinkdim(reference), _mean_chi(reference),
                        reference_error...,
                    ))
                    for cutoff in cutoffs, maxdim in maxdims
                        println("L_side=$side_level cutoff=$cutoff maxdim=$maxdim: compressing reference")
                        samples = _timed_compression(reference, cutoff, maxdim, repetitions)
                        representative = samples[end].field
                        direct_error = _probe_diagnostics(representative, sys)
                        _row(summary, (
                            side_level, total_bits, N, "truncate", cutoff, maxdim, repetitions,
                            _median([sample.time_s for sample in samples]), minimum(sample.time_s for sample in samples),
                            _median([sample.bytes for sample in samples]), _median([sample.gc_time_s for sample in samples]),
                            maxlinkdim(sys.ρ), _mean_chi(sys.ρ), maxlinkdim(representative), _mean_chi(representative),
                            direct_error...,
                        ))
                    end
                    flush(summary)
                catch error
                    _row(errors, (side_level, total_bits, "case", sprint(showerror, error, catch_backtrace())))
                    flush(errors)
                    println(stderr, "L_side=$side_level failed; details recorded in errors.csv")
                end
            end
        end
    end
    println("Compression-sweep artefacts written to $output")
end

function print_usage()
    println("Usage: julia --project=. benchmark/extraction_scaling_2d/adjacency_compression_sweep.jl --output DIRECTORY [options]")
    println("  --side-levels 4,6,8,10")
    println("  --cutoffs 1e-14,1e-12,1e-10,1e-8")
    println("  --maxdims 16,32,64,128,256")
    println("  --repetitions N")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_arguments(ARGS)
    isnothing(configuration) ? print_usage() : run_sweep(; configuration...)
end
