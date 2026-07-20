#!/usr/bin/env julia

using Dates
using ITensors, ITensorMPS
using LinearAlgebra
using MPO_MeanField
using Printf
using SHA

const DEFAULT_LEVELS = (4, 6, 8, 10, 12, 14)
const DEFAULT_SOURCES = (:smooth, :sp2_gapped)

function _csv_escape(value)
    string_value = string(value)
    return string('"', replace(string_value, '"' => "\"\""), '"')
end

function _write_csv_row(io, values)
    println(io, join(_csv_escape.(values), ','))
end

function _median(values::Vector{<:Real})
    ordered = sort(values)
    return ordered[cld(length(ordered), 2)]
end

function _bond_dimensions(mpo::MPO)
    length(mpo) <= 1 && return Int[]
    return [dim(commonind(mpo[site], mpo[site + 1])) for site in 1:(length(mpo) - 1)]
end

function _mean_bond_dimension(mpo::MPO)
    dimensions = _bond_dimensions(mpo)
    return isempty(dimensions) ? 1.0 : sum(dimensions) / length(dimensions)
end

function _dense_matrix(mpo::MPO, sys::System)
    N = 2^sys.params.L
    return [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states)
            for i in 1:N, j in 1:N]
end

function _mpo_relative_difference(left::MPO, right::MPO, params)
    difference = +(left, -right;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    numerator = sqrt(max(0.0, real(inner(difference, difference))))
    denominator = max(
        sqrt(max(0.0, real(inner(right, right)))),
        sqrt(eps(Float64)),
    )
    return numerator / denominator
end

function _hermiticity_residual(field::MPO, params)
    return _mpo_relative_difference(field, ITensors.dag(field), params)
end

function _sample_coefficients(field::MPO, sys::System, field_kind::Symbol)
    N = 2^sys.params.L
    positions = unique((1, max(1, div(N, 2)), N - 1))
    values = ComplexF64[]
    for position in positions
        value = if field_kind == :hartree
            MatrixChecker(field, sys.sites, position, position, sys.bra_states, sys.ket_states)
        elseif field_kind == :fock
            MatrixChecker(field, sys.sites, position, position + 1,
                sys.bra_states, sys.ket_states)
        else
            throw(ArgumentError("unknown field kind $field_kind"))
        end
        push!(values, value)
    end
    return join((@sprintf("%d:%.16e%+.16ei", position, real(value), imag(value))
                 for (position, value) in zip(positions, values)), ';')
end

function _smooth_density_system(L::Int)
    N = 2^L
    params = Parameters1D(
        L=L,
        t=-0.7,
        U=0.3,
        W=x -> 0.17 * cospi(Int(x) / (N - 1)) + 0.05 * sinpi(3Int(x) / (N - 1)),
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=64,
        density=0.5,
        purification_steps=45,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=5,
    )
    sys = System(params)
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        x -> 0.45 + 0.12 * cospi(Int(x) / (N - 1)),
        Float64, sys.sites, params.tci_tol,
    )
    bond = MPO_MeanField.diagonal_mpo_from_function(
        x -> Int(x) < N - 1 ? 0.04 * sinpi(3Int(x) / (N - 1)) : 0.0,
        Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L = sys.translations
    right = apply(bond, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    left = apply(T_L, ITensors.dag(bond); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.ρ = +(diagonal, right, left;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return sys, "deterministic_smooth_synthetic", "not_applicable"
end

function _sp2_gapped_density_system(L::Int)
    params = Parameters1D(
        L=L,
        t=-0.7,
        U=0.3,
        W=x -> iseven(Int(x)) ? 0.6 : -0.6,
        S=nothing,
        tci_tol=1e-10,
        itensors_tol=1e-12,
        itensors_maxdim=128,
        density=0.5,
        purification_steps=50,
        scf_mixing=0.5,
        scf_tol=0.1,
        scf_max_iterations=5,
    )
    sys = System(params)
    spectral_bounds = (-2.5, 2.5)
    rho0 = construct_rho_0(sys, params, spectral_bounds...; method=:sp2)
    result = perform_purification(
        rho0,
        params;
        method=:sp2,
        verbose=0,
        spectral_bounds=spectral_bounds,
    )
    result.converged || error(
        "SP2 preparation did not converge: $(result.termination_reason), " *
        "trace_error=$(result.trace_error), idempotency=$(result.idempotency_residual)",
    )
    sys.ρ = result.rho
    details = @sprintf(
        "sp2:iterations=%d;trace_error=%.3e;idempotency=%.3e;max_chi=%d",
        result.iterations,
        something(result.trace_error, NaN),
        result.idempotency_residual,
        result.final_bond_dimension,
    )
    return sys, "noninteracting_gapped_sp2", details
end

function _density_system(source::Symbol, L::Int)
    source == :smooth && return _smooth_density_system(L)
    source == :sp2_gapped && return _sp2_gapped_density_system(L)
    throw(ArgumentError("unknown density source $source"))
end

function _extractor(field_kind::Symbol, method::Symbol)
    if field_kind == :hartree
        method == :tci && return extract_hartree_mpo_1d
        method == :binary_carry && return extract_hartree_mpo_binary_carry_1d
    elseif field_kind == :fock
        method == :tci && return extract_fock_mpo_1d
        method == :binary_carry && return extract_fock_mpo_binary_carry_1d
    end
    throw(ArgumentError("unsupported field/method combination $field_kind/$method"))
end

function _benchmark_extractor(extractor::Function, sys::System, warmups::Int, repetitions::Int)
    warmups >= 0 || throw(ArgumentError("warmups must be nonnegative"))
    repetitions > 0 || throw(ArgumentError("repetitions must be positive"))
    for _ in 1:warmups
        extractor(sys)
    end
    timings = NamedTuple[]
    for repetition in 1:repetitions
        GC.gc(true)
        timing = @timed extractor(sys)
        field = timing.value
        push!(timings, (
            repetition=repetition,
            time_s=timing.time,
            allocations_bytes=timing.bytes,
            gc_time_s=timing.gctime,
            field_max_bond_dimension=maxlinkdim(field),
            field_mean_bond_dimension=_mean_bond_dimension(field),
            field_hermiticity_residual=_hermiticity_residual(field, sys.params),
            field=field,
        ))
    end
    return timings
end

function _write_metadata(output_dir::String, levels, sources, warmups, repetitions)
    project = Base.active_project()
    package_root = normpath(joinpath(@__DIR__, ".."))
    # Cluster source trees are commonly deployed with rsync and without `.git`.
    # Do not invoke Git (or emit its stderr) when provenance is unavailable.
    package_commit = "unavailable"
    if isdir(joinpath(package_root, ".git"))
        try
            package_commit = readchomp(pipeline(
                `git -C $package_root rev-parse HEAD`, stderr=devnull,
            ))
        catch
            # Git metadata is optional for a deployed cluster source tree.
        end
    end
    project_hash = isfile(project) ? bytes2hex(sha256(read(project))) : "unavailable"
    manifest_path = joinpath(dirname(project), "Manifest.toml")
    manifest_hash = isfile(manifest_path) ? bytes2hex(sha256(read(manifest_path))) : "unavailable"
    open(joinpath(output_dir, "metadata.txt"), "w") do io
        println(io, "started_at_utc=", Dates.now(Dates.UTC))
        println(io, "package_root=", package_root)
        println(io, "package_git_commit=", package_commit)
        println(io, "active_project=", project)
        println(io, "project_sha256=", project_hash)
        println(io, "manifest_path=", manifest_path)
        println(io, "manifest_sha256=", manifest_hash)
        println(io, "julia_version=", VERSION)
        println(io, "cpu_name=", Sys.CPU_NAME)
        println(io, "cpu_threads=", Sys.CPU_THREADS)
        println(io, "julia_threads=", Threads.nthreads())
        println(io, "blas_threads=", BLAS.get_num_threads())
        println(io, "total_memory_bytes=", Sys.total_memory())
        println(io, "kernel=", Sys.KERNEL)
        println(io, "machine=", Sys.MACHINE)
        println(io, "levels=", join(levels, ','))
        println(io, "sources=", join(string.(sources), ','))
        println(io, "warmups_per_extractor=", warmups)
        println(io, "timed_repetitions_per_extractor=", repetitions)
        println(io, "peak_rss_note=process-global RSS is intentionally not reported per method; use the wrapper's OS-level time output for the one-process peak")
    end
end

function _write_readme(output_dir::String)
    open(joinpath(output_dir, "README.txt"), "w") do io
        println(io, "Extraction benchmark artefacts")
        println(io, "")
        println(io, "samples.csv: one timed extraction per row; allocation values are cumulative allocations from @timed, not peak memory.")
        println(io, "summary.csv: median/minimum time and equivalence diagnostics per source/L/field/method.")
        println(io, "errors.csv: preparation or measurement failures; the remaining cases continue.")
        println(io, "metadata.txt: software, hardware, thread, truncation, and run configuration details.")
        println(io, "")
        println(io, "For L <= 4, dense_error_opnorm compares the two field matrices exactly.")
        println(io, "For larger L, relative_mpo_difference, sampled coefficients, and Hermiticity are MPO-level diagnostics; no large MPO is densified.")
        println(io, "The OS wrapper writes process_time.txt. Its peak RSS covers the whole single Julia process, not an individual extractor.")
    end
end

function _parse_symbols(value::String)
    return Tuple(Symbol(strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end

function _parse_levels(value::String)
    return Tuple(parse(Int, strip(item)) for item in split(value, ',') if !isempty(strip(item)))
end

function parse_arguments(arguments)
    output_dir = nothing
    levels = DEFAULT_LEVELS
    sources = DEFAULT_SOURCES
    warmups = 1
    repetitions_small = 5
    repetitions_large = 3
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (output_dir = arguments[index + 1]; index += 2; continue)
        argument == "--levels" && (levels = _parse_levels(arguments[index + 1]); index += 2; continue)
        argument == "--sources" && (sources = _parse_symbols(arguments[index + 1]); index += 2; continue)
        argument == "--warmups" && (warmups = parse(Int, arguments[index + 1]); index += 2; continue)
        argument == "--repetitions-small" && (repetitions_small = parse(Int, arguments[index + 1]); index += 2; continue)
        argument == "--repetitions-large" && (repetitions_large = parse(Int, arguments[index + 1]); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument; use --help"))
    end
    isnothing(output_dir) && throw(ArgumentError("--output DIRECTORY is required"))
    all(level -> level >= 1, levels) || throw(ArgumentError("all levels must be positive"))
    warmups >= 0 || throw(ArgumentError("--warmups must be nonnegative"))
    repetitions_small > 0 || throw(ArgumentError("--repetitions-small must be positive"))
    repetitions_large > 0 || throw(ArgumentError("--repetitions-large must be positive"))
    return (
        output_dir=abspath(output_dir),
        levels=levels,
        sources=sources,
        warmups=warmups,
        repetitions_small=repetitions_small,
        repetitions_large=repetitions_large,
    )
end

function print_usage(io::IO=stdout)
    println(io, "Usage: julia --project=. benchmark/extraction_scaling.jl --output DIRECTORY [options]")
    println(io, "  --levels 4,6,8,10,12,14")
    println(io, "  --sources smooth,sp2_gapped")
    println(io, "  --warmups N")
    println(io, "  --repetitions-small N   (levels <= 10; default 5)")
    println(io, "  --repetitions-large N   (levels >= 12; default 3)")
end

"""
    run_extraction_scaling_benchmark(; output_dir, ...)

Run the complete extraction comparison in one Julia process and preserve raw
samples plus analysis artefacts. Density preparation is outside extractor
timings. This is a performance/equivalence benchmark, not an SCF or physical
large-system calculation.
"""
function run_extraction_scaling_benchmark(
    ; output_dir::String,
    levels=DEFAULT_LEVELS,
    sources=DEFAULT_SOURCES,
    warmups::Int=1,
    repetitions_small::Int=5,
    repetitions_large::Int=3,
)
    mkpath(output_dir)
    _write_metadata(output_dir, levels, sources, warmups, repetitions_small)
    _write_readme(output_dir)
    sample_path = joinpath(output_dir, "samples.csv")
    summary_path = joinpath(output_dir, "summary.csv")
    error_path = joinpath(output_dir, "errors.csv")
    open(sample_path, "w") do sample_io
        _write_csv_row(sample_io, (
            "source", "L", "N", "field", "method", "repetition", "time_s",
            "allocations_bytes", "gc_time_s", "rho_max_chi", "rho_mean_chi",
            "field_max_chi", "field_mean_chi", "field_hermiticity_residual",
        ))
        open(summary_path, "w") do summary_io
            _write_csv_row(summary_io, (
                "source", "source_details", "L", "N", "field", "method",
                "repetitions", "median_time_s", "minimum_time_s",
                "median_allocations_bytes", "median_gc_time_s", "rho_max_chi",
                "rho_mean_chi", "median_field_max_chi", "median_field_mean_chi",
                "relative_mpo_difference_to_tci", "dense_error_opnorm",
                "tci_sample_coefficients", "carry_sample_coefficients",
                "tci_hermiticity_residual", "carry_hermiticity_residual",
            ))
            open(error_path, "w") do error_io
                _write_csv_row(error_io, ("source", "L", "stage", "error"))
                for source in sources, L in levels
                    N = 2^L
                    repetitions = L <= 10 ? repetitions_small : repetitions_large
                    println("source=$(source) L=$(L) N=$(N): preparing density")
                    sys = nothing
                    source_label = ""
                    source_details = ""
                    try
                        sys, source_label, source_details = _density_system(source, L)
                        fields = Dict{Tuple{Symbol,Symbol},MPO}()
                        timings = Dict{Tuple{Symbol,Symbol},Vector{NamedTuple}}()
                        for field_kind in (:hartree, :fock), method in (:tci, :binary_carry)
                            println("source=$(source) L=$(L) field=$(field_kind) method=$(method): timing $(repetitions) repetitions")
                            measurements = _benchmark_extractor(
                                _extractor(field_kind, method), sys, warmups, repetitions,
                            )
                            timings[(field_kind, method)] = measurements
                            fields[(field_kind, method)] = measurements[end].field
                            for measurement in measurements
                                _write_csv_row(sample_io, (
                                    source_label, L, N, field_kind, method,
                                    measurement.repetition, measurement.time_s,
                                    measurement.allocations_bytes, measurement.gc_time_s,
                                    maxlinkdim(sys.ρ), _mean_bond_dimension(sys.ρ),
                                    measurement.field_max_bond_dimension,
                                    measurement.field_mean_bond_dimension,
                                    measurement.field_hermiticity_residual,
                                ))
                            end
                        end
                        for field_kind in (:hartree, :fock)
                            tci = fields[(field_kind, :tci)]
                            carry = fields[(field_kind, :binary_carry)]
                            relative_error = _mpo_relative_difference(carry, tci, sys.params)
                            dense_error = L <= 4 ? opnorm(_dense_matrix(carry, sys) - _dense_matrix(tci, sys)) : NaN
                            tci_samples = _sample_coefficients(tci, sys, field_kind)
                            carry_samples = _sample_coefficients(carry, sys, field_kind)
                            for method in (:tci, :binary_carry)
                                measurements = timings[(field_kind, method)]
                                _write_csv_row(summary_io, (
                                    source_label, source_details, L, N, field_kind, method,
                                    repetitions,
                                    _median([item.time_s for item in measurements]),
                                    minimum(item.time_s for item in measurements),
                                    _median([item.allocations_bytes for item in measurements]),
                                    _median([item.gc_time_s for item in measurements]),
                                    maxlinkdim(sys.ρ), _mean_bond_dimension(sys.ρ),
                                    _median([item.field_max_bond_dimension for item in measurements]),
                                    _median([item.field_mean_bond_dimension for item in measurements]),
                                    relative_error, dense_error, tci_samples, carry_samples,
                                    _median([item.field_hermiticity_residual for item in timings[(field_kind, :tci)]]),
                                    _median([item.field_hermiticity_residual for item in timings[(field_kind, :binary_carry)]]),
                                ))
                            end
                        end
                        flush(sample_io)
                        flush(summary_io)
                    catch error
                        _write_csv_row(error_io, (source, L, "case", sprint(showerror, error, catch_backtrace())))
                        flush(error_io)
                        println(stderr, "source=$(source) L=$(L) failed; details recorded in errors.csv")
                    end
                end
            end
        end
    end
    println("Benchmark artefacts written to $output_dir")
    return output_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    configuration = parse_arguments(ARGS)
    if isnothing(configuration)
        print_usage()
    else
        run_extraction_scaling_benchmark(; configuration...)
    end
end
