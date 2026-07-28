#!/usr/bin/env julia

using Dates
using Printf
include("common.jl")

csv(value) = "\"" * replace(string(value), '"' => "\"\"") * "\""
row(io, values) = println(io, join(csv.(values), ','))

function parse_arguments(arguments)
    config = (
        output=nothing,
        side_level=10,
        cutoff=1e-10,
        maxdim=256,
        steps=50,
    )
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument == "--output" && (config=merge(config, (output=arguments[index + 1],)); index += 2; continue)
        argument == "--side-level" && (config=merge(config, (side_level=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--cutoff" && (config=merge(config, (cutoff=parse(Float64, arguments[index + 1]),)); index += 2; continue)
        argument == "--maxdim" && (config=merge(config, (maxdim=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--steps" && (config=merge(config, (steps=parse(Int, arguments[index + 1]),)); index += 2; continue)
        argument == "--help" && return nothing
        throw(ArgumentError("unknown argument $argument"))
    end
    isnothing(config.output) && throw(ArgumentError("--output is required"))
    1 <= config.side_level <= 10 || throw(ArgumentError("side level must lie in 1:10"))
    config.cutoff > 0 && config.maxdim > 0 && config.steps > 0 ||
        throw(ArgumentError("cutoff, maxdim, and steps must be positive"))
    config
end

function timed_field(operation)
    GC.gc(true)
    timing = @timed operation()
    field = timing.value
    return (
        field=field,
        time_s=timing.time,
        bytes=timing.bytes,
        gc_time_s=timing.gctime,
        max_chi=maxlinkdim(field),
        mean_chi=mean_bond_dimension(field),
    )
end

function run_profile(config)
    output = abspath(config.output)
    mkpath(output)
    preparation = @timed fixed_profile_system(
        config.side_level;
        cutoff=config.cutoff,
        maxdim=config.maxdim,
        steps=config.steps,
    )
    sys, purification, bounds = preparation.value

    println("Warming all mean-field component kernels outside timing")
    warm_hartree = extract_hartree_mpo_binary_carry_square_adjacency(
        sys; cutoff=0.0, maxdim=typemax(Int),
    )
    warm_hartree = nothing
    GC.gc(true)
    warm_horizontal = extract_fock_mpo_square_horizontal(sys)
    warm_vertical = extract_fock_mpo_square_vertical(sys)
    +(warm_horizontal, warm_vertical; cutoff=config.cutoff, maxdim=config.maxdim)
    warm_horizontal = warm_vertical = nothing
    GC.gc(true)
    warm_horizontal_carry = extract_fock_mpo_binary_carry_square_horizontal(sys)
    warm_vertical_carry = extract_fock_mpo_binary_carry_square_vertical(sys)
    +(warm_horizontal_carry, warm_vertical_carry;
        cutoff=config.cutoff, maxdim=config.maxdim)
    warm_horizontal_carry = warm_vertical_carry = nothing
    GC.gc(true)

    hartree = timed_field(() -> extract_hartree_mpo_binary_carry_square_adjacency(
        sys; cutoff=0.0, maxdim=typemax(Int),
    ))
    horizontal = timed_field(() -> extract_fock_mpo_square_horizontal(sys))
    vertical = timed_field(() -> extract_fock_mpo_square_vertical(sys))
    fock_sum = timed_field(() -> +(
        horizontal.field,
        vertical.field;
        cutoff=config.cutoff,
        maxdim=config.maxdim,
    ))
    horizontal_carry = timed_field(
        () -> extract_fock_mpo_binary_carry_square_horizontal(sys),
    )
    vertical_carry = timed_field(
        () -> extract_fock_mpo_binary_carry_square_vertical(sys),
    )
    carry_sum = timed_field(() -> +(
        horizontal_carry.field,
        vertical_carry.field;
        cutoff=config.cutoff,
        maxdim=config.maxdim,
    ))

    open(joinpath(output, "components.csv"), "w") do io
        row(io, (
            "component", "time_s", "allocations_bytes", "gc_time_s",
            "max_chi", "mean_chi", "relative_error_to_current",
        ))
        for (label, measurement, error) in (
            ("hartree_adjacency_untruncated", hartree, 0.0),
            ("fock_horizontal_tci", horizontal, 0.0),
            ("fock_vertical_tci", vertical, 0.0),
            ("fock_sum_tci", fock_sum, 0.0),
            (
                "fock_horizontal_binary_carry",
                horizontal_carry,
                relative_mpo_error(horizontal_carry.field, horizontal.field),
            ),
            (
                "fock_vertical_binary_carry",
                vertical_carry,
                relative_mpo_error(vertical_carry.field, vertical.field),
            ),
            (
                "fock_sum_binary_carry",
                carry_sum,
                relative_mpo_error(carry_sum.field, fock_sum.field),
            ),
        )
            row(io, (
                label, measurement.time_s, measurement.bytes,
                measurement.gc_time_s, measurement.max_chi,
                measurement.mean_chi, error,
            ))
        end
    end

    open(joinpath(output, "metadata.toml"), "w") do io
        println(io, "created_at = \"$(Dates.now(Dates.UTC))\"")
        println(io, "side_level = $(config.side_level)")
        println(io, "L_total = $(2config.side_level)")
        println(io, "N = $(2^(2config.side_level))")
        println(io, "cutoff = $(config.cutoff)")
        println(io, "maxdim = $(config.maxdim)")
        println(io, "preparation_time_s = $(preparation.time)")
        println(io, "preparation_allocations_bytes = $(preparation.bytes)")
        println(io, "purification_iterations = $(purification.iterations)")
        println(io, "purification_trace_error = $(purification.trace_error)")
        println(io, "purification_idempotency = $(purification.idempotency_residual)")
        println(io, "rho_max_chi = $(maxlinkdim(sys.ρ))")
        println(io, "spectral_lower = $(bounds[1])")
        println(io, "spectral_upper = $(bounds[2])")
    end
    println("Component profile written to $output")
end

if abspath(PROGRAM_FILE) == @__FILE__
    config = parse_arguments(ARGS)
    if isnothing(config)
        println("Usage: component_profile_2d.jl --output DIR [--side-level 10 --cutoff 1e-10 --maxdim 256]")
    else
        run_profile(config)
    end
end
