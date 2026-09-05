#!/usr/bin/env julia

"""Measure Hadamard probing error at selected entries without diagonalization.

The deterministic reference is the same Jackson--Chebyshev polynomial applied
to basis vectors for a fixed grid of sites and their right/up neighbours.
Consequently, differences measure probing error only: there is no dense exact
diagonalization, full projector, or complete Hadamard basis.
"""

using Dates
using LinearAlgebra
using MPO_MeanField
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) == 8 || error(
    "usage: diagnose_square_kpm_selected.jl CAMPAIGN_FILE TASK_INDEX " *
    "MOMENTS PROBES SEED SAMPLE_GRID OUTPUT_DIRECTORY BACKEND",
)

campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
moments = parse(Int, ARGS[3])
probes = parse(Int, ARGS[4])
seed = parse(Int, ARGS[5])
sample_grid = parse(Int, ARGS[6])
output = abspath(ARGS[7])
backend = Symbol(lowercase(ARGS[8]))

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
moments >= 2 || error("MOMENTS must be at least 2")
probes > 0 || error("PROBES must be positive")
sample_grid >= 2 || error("SAMPLE_GRID must be at least 2")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
1 <= task_index <= length(campaign.runs) || error("TASK_INDEX is outside campaign")
spec = campaign.runs[task_index]
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

function direct_initial_hamiltonian(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    rows = Int[]
    columns = Int[]
    values = Float64[]
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) :
        Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) :
        Float64(params.t[2](x, y))
    for x in 0:(side - 1), y in 0:(side - 1)
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
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, params.L)
            hopping = ty(x, y)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour, site))
            append!(values, (hopping, hopping))
        end
    end
    return sparse(rows, columns, values, N, N)
end

function selected_geometry(params::ParametersSquare, grid::Int)
    side = 2^div(params.L, 2)
    coordinates = unique(round.(Int, range(0, side - 1; length=grid)))
    sites = sort(unique(
        square_lattice_index(x, y, params.L)
        for x in coordinates for y in coordinates
    ))
    bonds = Tuple{Int,Int,Symbol}[]
    for site in sites
        neighbours = square_neighbours(site, params.L)
        !isnothing(neighbours.right) &&
            push!(bonds, (site, neighbours.right, :horizontal))
        !isnothing(neighbours.up) &&
            push!(bonds, (site, neighbours.up, :vertical))
    end
    columns = sort(unique(vcat(
        sites,
        [site for bond in bonds for site in (bond[1], bond[2])],
    )))
    return (; sites, bonds, columns)
end

N = 2^params.L
probes <= N || error("PROBES=$probes exceeds N=$N")
Ne = round(Int, params.density * N)
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
center = (lower + upper) / 2
halfwidth = (upper - lower) / 2
scaled_mu = -center / halfwidth
coefficients = projector_coefficients(moments, scaled_mu)
H = direct_initial_hamiltonian(params)
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, N, N)) / halfwidth
geometry = selected_geometry(params, sample_grid)
column_position = Dict(site => index for (index, site) in enumerate(geometry.columns))
row_position = Dict(site => index for (index, site) in enumerate(geometry.columns))

host_probes = probing_matrix(N, probes, :hadamard, seed)
basis = zeros(Float64, N, length(geometry.columns))
for (column, site) in enumerate(geometry.columns)
    basis[site, column] = 1.0
end
host_vectors = hcat(host_probes, basis)

device_name = "CPU"
device_total_memory = 0
device_free_memory_before = 0
transfer_time = 0.0
synchronize = () -> nothing
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    device_total_memory = CUDA.total_memory()
    device_free_memory_before = CUDA.free_memory()
    transfer = @timed begin
        H_backend = CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H)
        vectors_backend = CUDA.CuArray(host_vectors)
        (; H_backend, vectors_backend)
    end
    transfer_time = transfer.time
    H_backend = transfer.value.H_backend
    vectors_backend = transfer.value.vectors_backend
    synchronize = CUDA.synchronize
else
    H_backend = scaled_H
    vectors_backend = host_vectors
end

synchronize()
kpm_calculation = @timed begin
    filtered_backend = kpm_apply(
        H_backend,
        vectors_backend,
        coefficients;
        synchronize=synchronize,
    )
    synchronize()
    filtered_backend
end
filtered_backend = kpm_calculation.value

selected_filtered = backend == :cuda ?
    Array(filtered_backend[geometry.columns, :]) :
    filtered_backend[geometry.columns, :]
selected_probes = @view host_probes[geometry.columns, :]
probe_filtered = @view selected_filtered[:, 1:probes]
reference_filtered =
    @view selected_filtered[:, (probes + 1):size(selected_filtered, 2)]

density_estimate = Vector{Float64}(undef, length(geometry.sites))
density_reference = similar(density_estimate)
for (index, site) in enumerate(geometry.sites)
    row = row_position[site]
    density_estimate[index] =
        dot(@view(probe_filtered[row, :]), @view(selected_probes[row, :])) /
        probes
    density_reference[index] =
        reference_filtered[row, column_position[site]]
end

bond_estimate = Vector{Float64}(undef, length(geometry.bonds))
bond_reference = similar(bond_estimate)
for (index, (site, neighbour, _)) in enumerate(geometry.bonds)
    row_site = row_position[site]
    row_neighbour = row_position[neighbour]
    bond_estimate[index] = (
        dot(
            @view(probe_filtered[row_site, :]),
            @view(selected_probes[row_neighbour, :]),
        ) +
        dot(
            @view(probe_filtered[row_neighbour, :]),
            @view(selected_probes[row_site, :]),
        )
    ) / (2probes)
    bond_reference[index] = (
        reference_filtered[row_site, column_position[neighbour]] +
        reference_filtered[row_neighbour, column_position[site]]
    ) / 2
end

density_error = density_estimate - density_reference
bond_error = bond_estimate - bond_reference
horizontal = findall(bond -> bond[3] == :horizontal, geometry.bonds)
vertical = findall(bond -> bond[3] == :vertical, geometry.bonds)

# The full trace estimate is a Frobenius inner product. It does not require
# transferring the N-by-R filtered block back to the host.
trace_estimate = if backend == :cuda
    Float64(dot(
        vec(@view(filtered_backend[:, 1:probes])),
        vec(@view(vectors_backend[:, 1:probes])),
    )) / probes
else
    dot(
        vec(@view(filtered_backend[:, 1:probes])),
        vec(@view(vectors_backend[:, 1:probes])),
    ) / probes
end

summary = Dict(
    "campaign" => campaign.name,
    "label" => spec.label,
    "diagnostic" => "selected_entry_hadamard_kpm",
    "backend" => string(backend),
    "device_name" => device_name,
    "device_total_memory_bytes" => device_total_memory,
    "device_free_memory_before_bytes" => device_free_memory_before,
    "matrix_dimension" => N,
    "target_particles" => Ne,
    "moments" => moments,
    "probes" => probes,
    "seed" => seed,
    "sample_grid" => sample_grid,
    "sampled_sites" => length(geometry.sites),
    "sampled_bonds" => length(geometry.bonds),
    "reference_basis_vectors" => length(geometry.columns),
    "spectral_lower" => lower,
    "spectral_upper" => upper,
    "trace" => trace_estimate,
    "trace_error" => abs(trace_estimate - Ne),
    "selected_density_max_abs_error" => maximum(abs, density_error),
    "selected_density_rms_error" =>
        norm(density_error) / sqrt(length(density_error)),
    "selected_horizontal_bond_max_abs_error" =>
        maximum(abs, @view bond_error[horizontal]),
    "selected_horizontal_bond_rms_error" =>
        norm(@view bond_error[horizontal]) / sqrt(length(horizontal)),
    "selected_vertical_bond_max_abs_error" =>
        maximum(abs, @view bond_error[vertical]),
    "selected_vertical_bond_rms_error" =>
        norm(@view bond_error[vertical]) / sqrt(length(vertical)),
    "transfer_to_device_time_s" => transfer_time,
    "kpm_time_s" => kpm_calculation.time,
    "kpm_allocations_bytes" => kpm_calculation.bytes,
    "finished_at" => string(now(UTC)),
)

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, summary)
end
open(joinpath(output, "selected_density.csv"), "w") do io
    write_csv_row(io, ("site", "x", "y", "reference", "estimate", "error"))
    for (index, site) in enumerate(geometry.sites)
        x, y = square_lattice_decoder(site - 1, params.L)
        write_csv_row(io, (
            site, x, y, density_reference[index], density_estimate[index],
            density_error[index],
        ))
    end
end
open(joinpath(output, "selected_bonds.csv"), "w") do io
    write_csv_row(io, (
        "site", "neighbour", "orientation", "reference", "estimate", "error",
    ))
    for (index, bond) in enumerate(geometry.bonds)
        write_csv_row(io, (
            bond[1], bond[2], bond[3], bond_reference[index],
            bond_estimate[index], bond_error[index],
        ))
    end
end

println("Selected-entry KPM diagnostic complete: $output")
println("  N=$N moments=$moments probes=$probes")
println("  sampled sites=$(length(geometry.sites)) bonds=$(length(geometry.bonds))")
@printf("  selected density max/RMS error = %.6e / %.6e\n",
    summary["selected_density_max_abs_error"],
    summary["selected_density_rms_error"])

