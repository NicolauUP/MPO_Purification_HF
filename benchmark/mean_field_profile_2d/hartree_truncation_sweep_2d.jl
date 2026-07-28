#!/usr/bin/env julia

using Dates
include("common.jl")

csv(value) = "\"" * replace(string(value), '"' => "\"\"") * "\""
row(io, values) = println(io, join(csv.(values), ','))
parse_list(value, type) = Tuple(parse(type, strip(item)) for item in split(value, ','))

function parse_arguments(arguments)
    config = (
        output=nothing,
        side_level=10,
        preparation_cutoff=1e-10,
        preparation_maxdim=256,
        cutoffs=(1e-12, 1e-10, 1e-8),
        maxdims=(64, 128, 256, 384, 512),
        steps=50,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (config=merge(config, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (config=merge(config, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--preparation-cutoff" && (config=merge(config, (preparation_cutoff=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--preparation-maxdim" && (config=merge(config, (preparation_maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--cutoffs" && (config=merge(config, (cutoffs=parse_list(arguments[index + 1], Float64),)); index += 2; continue)
        argument == "--maxdims" && (config=merge(config, (maxdims=parse_list(arguments[index + 1], Int),)); index += 2; continue)
        argument == "--steps" && (config=merge(config, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument"))
    end
    isnothing(config.output) && throw(ArgumentError("--output is required"))
    1 <= config.side_level <= 10 || throw(ArgumentError("side level must lie in 1:10"))
    all(>(0), config.cutoffs) && all(>(0), config.maxdims) ||
        throw(ArgumentError("cutoffs and maxdims must be positive"))
    config
end

function measure_field(operation)
    GC.gc(true)
    timing = @timed operation()
    field = timing.value
    probe_max, probe_mean = direct_hartree_probe_error(field, CURRENT_SYSTEM[])
    return (
        field=field,
        time_s=timing.time,
        bytes=timing.bytes,
        gc_time_s=timing.gctime,
        max_chi=maxlinkdim(field),
        mean_chi=mean_bond_dimension(field),
        probe_max=probe_max,
        probe_mean=probe_mean,
    )
end

const CURRENT_SYSTEM = Ref{Any}(nothing)

function run_sweep(config)
    output = abspath(config.output)
    mkpath(output)
    preparation = @timed fixed_profile_system(
        config.side_level;
        cutoff=config.preparation_cutoff,
        maxdim=config.preparation_maxdim,
        steps=config.steps,
    )
    sys, purification, _ = preparation.value
    CURRENT_SYSTEM[] = sys

    println("Warming untruncated adjacency application outside timing")
    extract_hartree_mpo_binary_carry_square_adjacency(sys)
    reference = measure_field(
        () -> extract_hartree_mpo_binary_carry_square_adjacency(sys),
    )
    open(joinpath(output, "summary.csv"), "w") do io
        row(io, (
            "policy", "cutoff", "maxdim", "time_s", "allocations_bytes",
            "gc_time_s", "field_max_chi", "field_mean_chi",
            "relative_error_to_untruncated", "direct_probe_max_abs_error",
            "direct_probe_mean_abs_error",
        ))
        row(io, (
            "untruncated", 0.0, "unbounded", reference.time_s,
            reference.bytes, reference.gc_time_s, reference.max_chi,
            reference.mean_chi, 0.0, reference.probe_max,
            reference.probe_mean,
        ))
        for cutoff in config.cutoffs, maxdim in config.maxdims
            println("cutoff=$cutoff maxdim=$maxdim")
            truncated_adjacency_hartree(
                sys;
                cutoff=cutoff,
                maxdim=maxdim,
            )
            candidate = measure_field(
                () -> truncated_adjacency_hartree(
                    sys;
                    cutoff=cutoff,
                    maxdim=maxdim,
                ),
            )
            row(io, (
                "truncated_apply", cutoff, maxdim, candidate.time_s,
                candidate.bytes, candidate.gc_time_s, candidate.max_chi,
                candidate.mean_chi,
                relative_mpo_error(candidate.field, reference.field),
                candidate.probe_max, candidate.probe_mean,
            ))
            flush(io)
        end
    end
    open(joinpath(output, "metadata.toml"), "w") do io
        println(io, "created_at = \"$(Dates.now(Dates.UTC))\"")
        println(io, "side_level = $(config.side_level)")
        println(io, "L_total = $(2config.side_level)")
        println(io, "N = $(2^(2config.side_level))")
        println(io, "preparation_cutoff = $(config.preparation_cutoff)")
        println(io, "preparation_maxdim = $(config.preparation_maxdim)")
        println(io, "preparation_time_s = $(preparation.time)")
        println(io, "purification_iterations = $(purification.iterations)")
        println(io, "rho_max_chi = $(maxlinkdim(sys.ρ))")
        println(io, "cutoffs = $(collect(config.cutoffs))")
        println(io, "maxdims = $(collect(config.maxdims))")
    end
    println("Hartree truncation sweep written to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    config = parse_arguments(ARGS)
    if isnothing(config)
        println("Usage: hartree_truncation_sweep_2d.jl --output DIR [options]")
    else
        run_sweep(config)
    end
end
