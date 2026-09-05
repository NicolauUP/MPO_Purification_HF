#!/usr/bin/env julia

"""All-site nested-Hadamard KPM convergence without exact diagonalization.

The square-lattice matrix rows retain the project's interleaved quantics
ordering, y₀,x₀,y₁,x₁,... . Hadamard prefixes are therefore coordinate-aware.
The largest requested prefix is an operational reference for probing error;
all configurations use the same Jackson--Chebyshev polynomial.
"""

using Dates
using LinearAlgebra
using MPO_MeanField
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 2:3 || error(
    "usage: diagnose_square_kpm_all_nested.jl " *
    "CAMPAIGN_FILE OUTPUT_DIRECTORY [BACKEND]",
)
campaign_file = abspath(ARGS[1])
output = abspath(ARGS[2])
backend = length(ARGS) == 3 ? Symbol(lowercase(ARGS[3])) : :cpu

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
length(campaign.runs) == 1 || error("diagnostic expects one campaign run")
spec = only(campaign.runs)
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

const MOMENTS = 800
const PRIMARY_SEED = 510578
const REFERENCE_PROBES = 8192
const PROBE_CONFIGURATIONS = [
    (512, PRIMARY_SEED),
    (1024, PRIMARY_SEED),
    (2048, PRIMARY_SEED),
    (4096, PRIMARY_SEED),
    (8192, PRIMARY_SEED),
    (1024, 20260730),
    (2048, 20260730),
    (4096, 20260730),
]

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

function direct_initial_hamiltonian_and_bonds(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    rows = Int[]
    columns = Int[]
    values = Float64[]
    bonds = Tuple{Int,Int,Symbol}[]
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) :
        Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) :
        Float64(params.t[2](x, y))
    for y in 0:(side - 1), x in 0:(side - 1)
        site = square_lattice_index(x, y, params.L)
        onsite = (isnothing(params.W) ? 0.0 : Float64(params.W(x, y))) +
                 (isnothing(params.S) ? 0.0 : Float64(params.S(x, y)))
        push!(rows, site); push!(columns, site); push!(values, onsite)
        if x < side - 1
            neighbour = square_lattice_index(x + 1, y, params.L)
            hopping = tx(x, y)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour, site))
            append!(values, (hopping, hopping))
            push!(bonds, (site, neighbour, :horizontal))
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, params.L)
            hopping = ty(x, y)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour, site))
            append!(values, (hopping, hopping))
            push!(bonds, (site, neighbour, :vertical))
        end
    end
    return sparse(rows, columns, values, N, N), bonds
end

function local_estimates(filtered, probes, bonds)
    R = size(probes, 2)
    density = vec(mean(filtered .* probes; dims=2))
    bond_order = Vector{Float64}(undef, length(bonds))
    for (index, (site, neighbour, _)) in enumerate(bonds)
        bond_order[index] = (
            dot(@view(filtered[site, :]), @view(probes[neighbour, :])) +
            dot(@view(filtered[neighbour, :]), @view(probes[site, :]))
        ) / (2R)
    end
    return (; density, bond_order)
end

function hartree_field(density, params::ParametersSquare)
    field = Vector{Float64}(undef, length(density))
    for site in eachindex(density)
        neighbours = square_neighbours(site, params.L)
        field[site] = params.U * sum(
            density[neighbour]
            for neighbour in values(neighbours) if !isnothing(neighbour)
        )
    end
    return field
end

function checkerboard_order(density, params::ParametersSquare)
    side = 2^div(params.L, 2)
    return sum(0:(side - 1)) do y
        sum(0:(side - 1)) do x
            sign = iseven(x + y) ? 1.0 : -1.0
            sign * density[square_lattice_index(x, y, params.L)]
        end
    end / length(density)
end

function one_body_energy(H, density, bonds, bond_order)
    onsite = sum(H[site, site] * density[site] for site in eachindex(density))
    hopping = sum(
        2H[site, neighbour] * bond_order[index]
        for (index, (site, neighbour, _)) in enumerate(bonds)
    )
    return onsite + hopping
end

function error_metrics(values, reference)
    difference = values - reference
    return (
        maximum(abs, difference),
        norm(difference) / sqrt(length(difference)),
    )
end

N = 2^params.L
side = 2^div(params.L, 2)
N == 65_536 || error("this experiment expects N=65,536, got N=$N")
REFERENCE_PROBES <= N || error("reference probe count exceeds N")
Ne = round(Int, params.density * N)
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
center = (lower + upper) / 2
halfwidth = (upper - lower) / 2
scaled_mu = -center / halfwidth
coefficients = projector_coefficients(MOMENTS, scaled_mu)
H, bonds = direct_initial_hamiltonian_and_bonds(params)
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, N, N)) / halfwidth
horizontal = findall(bond -> bond[3] == :horizontal, bonds)
vertical = findall(bond -> bond[3] == :vertical, bonds)

device_name = "CPU"
device_total_memory = 0
device_free_memory_before = 0
synchronize = () -> nothing
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    device_total_memory = CUDA.total_memory()
    device_free_memory_before = CUDA.free_memory()
    H_backend = CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H)
    synchronize = CUDA.synchronize
else
    H_backend = scaled_H
end

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
metadata = Dict(
    "campaign" => campaign.name,
    "label" => spec.label,
    "diagnostic" => "all_site_nested_coordinate_hadamard_kpm",
    "backend" => string(backend),
    "device_name" => device_name,
    "device_total_memory_bytes" => device_total_memory,
    "device_free_memory_before_bytes" => device_free_memory_before,
    "side" => side,
    "matrix_dimension" => N,
    "target_particles" => Ne,
    "moments" => MOMENTS,
    "spectral_lower" => lower,
    "spectral_upper" => upper,
    "primary_seed" => PRIMARY_SEED,
    "operational_reference_probes" => REFERENCE_PROBES,
    "probe_ordering" => "project_square_interleaved_coordinate_bits",
    "reference_scope" => "probing convergence at fixed KPM polynomial",
    "probe_configurations" => [
        "R=$(configuration[1]),seed=$(configuration[2])"
        for configuration in PROBE_CONFIGURATIONS
    ],
    "created_at" => string(now(UTC)),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

raw_summary_path = joinpath(output, "raw_summary.csv")
open(raw_summary_path, "w") do io
    write_csv_row(io, (
        "probes", "seed", "kpm_time_s", "transfer_to_device_time_s",
        "transfer_to_host_time_s", "measurement_time_s", "trace",
        "trace_error", "checkerboard_order", "one_body_energy",
    ))
end

estimates = Dict{Tuple{Int,Int},NamedTuple}()
for (probes, seed) in PROBE_CONFIGURATIONS
    case_name = "m$(MOMENTS)_hadamard_r$(probes)_seed$(seed)"
    case_output = joinpath(output, case_name)
    mkpath(case_output)
    println("All-site coordinate-Hadamard KPM: R=$probes seed=$seed")
    flush(stdout)

    host_probes = probing_matrix(N, probes, :hadamard, seed)
    transfer_to_device = @timed begin
        probes_backend = backend == :cuda ?
            CUDA.CuArray(host_probes) : host_probes
        synchronize()
        probes_backend
    end
    probes_backend = transfer_to_device.value

    synchronize()
    kpm_calculation = @timed begin
        value = kpm_apply(
            H_backend, probes_backend, coefficients;
            synchronize=synchronize,
        )
        synchronize()
        value
    end
    filtered_backend = kpm_calculation.value
    transfer_to_host = @timed begin
        value = backend == :cuda ? Array(filtered_backend) : filtered_backend
        value
    end
    filtered = transfer_to_host.value
    measurement = @timed local_estimates(filtered, host_probes, bonds)
    density = measurement.value.density
    bond_order = measurement.value.bond_order
    hartree = hartree_field(density, params)
    fock = -params.U .* bond_order
    trace_value = sum(density)
    order = checkerboard_order(density, params)
    energy = one_body_energy(H, density, bonds, bond_order)
    estimates[(probes, seed)] = (;
        density, bond_order, hartree, fock, trace_value, order, energy,
    )

    open(raw_summary_path, "a") do io
        write_csv_row(io, (
            probes, seed, kpm_calculation.time, transfer_to_device.time,
            transfer_to_host.time, measurement.time, trace_value,
            abs(trace_value - Ne), order, energy,
        ))
    end
    open(joinpath(case_output, "density.csv"), "w") do io
        write_csv_row(io, ("site", "x", "y", "density", "hartree"))
        for y in 0:(side - 1), x in 0:(side - 1)
            site = square_lattice_index(x, y, params.L)
            write_csv_row(io, (site, x, y, density[site], hartree[site]))
        end
    end
    open(joinpath(case_output, "bonds.csv"), "w") do io
        write_csv_row(io, (
            "site", "neighbour", "orientation", "bond_order", "fock",
        ))
        for (index, bond) in enumerate(bonds)
            write_csv_row(io, (
                bond[1], bond[2], bond[3], bond_order[index], fock[index],
            ))
        end
    end

    @printf(
        "  Tr=%.12g checkerboard=%.6e KPM=%.3fs measurement=%.3fs\n",
        trace_value, order, kpm_calculation.time, measurement.time,
    )
    flush(stdout)

    host_probes = nothing
    probes_backend = nothing
    filtered_backend = nothing
    filtered = nothing
    GC.gc()
    backend == :cuda && CUDA.reclaim()
end

reference = estimates[(REFERENCE_PROBES, PRIMARY_SEED)]
convergence_path = joinpath(output, "convergence_to_r8192.csv")
open(convergence_path, "w") do io
    write_csv_row(io, (
        "probes", "seed", "trace_difference", "checkerboard_difference",
        "energy_difference", "density_max_abs_difference",
        "density_rms_difference", "hartree_max_abs_difference",
        "hartree_rms_difference", "horizontal_fock_max_abs_difference",
        "horizontal_fock_rms_difference",
        "vertical_fock_max_abs_difference",
        "vertical_fock_rms_difference",
    ))
    for (probes, seed) in PROBE_CONFIGURATIONS
        value = estimates[(probes, seed)]
        density_error = error_metrics(value.density, reference.density)
        hartree_error = error_metrics(value.hartree, reference.hartree)
        horizontal_fock_error = error_metrics(
            @view(value.fock[horizontal]), @view(reference.fock[horizontal]),
        )
        vertical_fock_error = error_metrics(
            @view(value.fock[vertical]), @view(reference.fock[vertical]),
        )
        write_csv_row(io, (
            probes, seed, value.trace_value - reference.trace_value,
            value.order - reference.order, value.energy - reference.energy,
            density_error..., hartree_error..., horizontal_fock_error...,
            vertical_fock_error...,
        ))
    end
end

seed_path = joinpath(output, "seed_differences.csv")
open(seed_path, "w") do io
    write_csv_row(io, (
        "probes", "seed_a", "seed_b", "trace_difference",
        "density_max_abs_difference", "density_rms_difference",
        "hartree_max_abs_difference", "hartree_rms_difference",
        "horizontal_fock_max_abs_difference",
        "horizontal_fock_rms_difference",
        "vertical_fock_max_abs_difference",
        "vertical_fock_rms_difference",
    ))
    for probes in (1024, 2048, 4096)
        a = estimates[(probes, PRIMARY_SEED)]
        b = estimates[(probes, 20260730)]
        density_error = error_metrics(a.density, b.density)
        hartree_error = error_metrics(a.hartree, b.hartree)
        horizontal_fock_error = error_metrics(
            @view(a.fock[horizontal]), @view(b.fock[horizontal]),
        )
        vertical_fock_error = error_metrics(
            @view(a.fock[vertical]), @view(b.fock[vertical]),
        )
        write_csv_row(io, (
            probes, PRIMARY_SEED, 20260730,
            a.trace_value - b.trace_value, density_error..., hartree_error...,
            horizontal_fock_error..., vertical_fock_error...,
        ))
    end
end

println("All-site nested-Hadamard comparison complete: $output")
println("  operational reference: R=$REFERENCE_PROBES seed=$PRIMARY_SEED")

