#!/usr/bin/env julia

"""Compare coordinate-interleaved Hadamard probes with a saved ED reference.

Rows remain in contiguous x-fastest physical storage order. Only the
Hadamard row code changes to an interleaving of the x and y coordinate bits.
The dense eigensystem is not recomputed; exact densities and bonds are loaded
from the completed 128x64 reference experiment.
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 2:3 || error(
    "usage: compare_rectangular_kpm_coordinate_hadamard.jl " *
    "REFERENCE_DIRECTORY OUTPUT_DIRECTORY [BACKEND]",
)
reference_directory = abspath(ARGS[1])
output = abspath(ARGS[2])
backend = length(ARGS) == 3 ? Symbol(lowercase(ARGS[3])) : :cpu
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

const NX = 128
const NY = 64
const N = NX * NY
const NE = div(N, 2)
const V2 = 0.5
const SEED_AMPLITUDE = 0.1
const TAU = sqrt(2.0) - 5.0 / 6.0
const SPECTRAL_LOWER = -6.35
const SPECTRAL_UPPER = 6.35
const MOMENTS = 800
const PROBE_CONFIGURATIONS = [
    (256, 510578),
    (512, 510578),
    (1024, 510578),
    (2048, 510578),
    (4096, 510578),
    (8192, 510578),
    (512, 20260730),
    (1024, 20260730),
    (2048, 20260730),
]

site_index(x, y) = x + NX * y + 1
tx(x) = -1.0 - V2 * cos(2π * TAU * (Float64(x) + 0.5))
ty(y) = -1.0 - V2 * cos(2π * TAU * (Float64(y) + 0.5))
seed_field(x, y) = iseven(x + y) ? SEED_AMPLITUDE : -SEED_AMPLITUDE

function coordinate_code(x::Int, y::Int)
    # Six paired low bits, followed by the unmatched highest x bit.
    code = 0
    for bit in 0:5
        code |= ((y >> bit) & 1) << (2bit)
        code |= ((x >> bit) & 1) << (2bit + 1)
    end
    code |= ((x >> 6) & 1) << 12
    return code
end

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
parse_csv_fields(line) = split(replace(chomp(line), "\"" => ""), ',')

function load_exact_reference(directory)
    metadata_path = joinpath(directory, "metadata.toml")
    density_path = joinpath(
        directory, "m800_hadamard_r4096_seed510578", "density.csv",
    )
    bonds_path = joinpath(
        directory, "m800_hadamard_r4096_seed510578", "bonds.csv",
    )
    isfile(metadata_path) || error("missing ED metadata: $metadata_path")
    isfile(density_path) || error("missing ED density: $density_path")
    isfile(bonds_path) || error("missing ED bonds: $bonds_path")
    metadata = TOML.parsefile(metadata_path)
    metadata["matrix_dimension"] == N ||
        error("reference matrix dimension does not equal $N")

    density = Vector{Float64}(undef, N)
    open(density_path) do io
        readline(io)
        count = 0
        for line in eachline(io)
            fields = parse_csv_fields(line)
            site = parse(Int, fields[1])
            density[site] = parse(Float64, fields[4])
            count += 1
        end
        count == N || error("expected $N density rows, found $count")
    end

    bonds = Tuple{Int,Int,Symbol}[]
    bond_order = Float64[]
    open(bonds_path) do io
        readline(io)
        for line in eachline(io)
            fields = parse_csv_fields(line)
            push!(bonds, (
                parse(Int, fields[1]),
                parse(Int, fields[2]),
                Symbol(fields[3]),
            ))
            push!(bond_order, parse(Float64, fields[4]))
        end
    end
    return (; metadata, density, bonds, bond_order)
end

function build_hamiltonian()
    rows = Int[]
    columns = Int[]
    values = Float64[]
    for y in 0:(NY - 1), x in 0:(NX - 1)
        site = site_index(x, y)
        push!(rows, site); push!(columns, site)
        push!(values, seed_field(x, y))
        if x < NX - 1
            neighbour = site_index(x + 1, y)
            hopping = tx(x)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour, site))
            append!(values, (hopping, hopping))
        end
        if y < NY - 1
            neighbour = site_index(x, y + 1)
            hopping = ty(y)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour, site))
            append!(values, (hopping, hopping))
        end
    end
    return sparse(rows, columns, values, N, N)
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

reference = load_exact_reference(reference_directory)
exact_density = reference.density
bonds = reference.bonds
exact_bond_order = reference.bond_order
horizontal = findall(bond -> bond[3] == :horizontal, bonds)
vertical = findall(bond -> bond[3] == :vertical, bonds)
exact_mu = (
    reference.metadata["exact_homo"] + reference.metadata["exact_lumo"]
) / 2

codes = Vector{Int}(undef, N)
for y in 0:(NY - 1), x in 0:(NX - 1)
    codes[site_index(x, y)] = coordinate_code(x, y)
end
sort(codes) == collect(0:(N - 1)) ||
    error("coordinate codes are not a permutation of 0:$(N - 1)")

center = (SPECTRAL_LOWER + SPECTRAL_UPPER) / 2
halfwidth = (SPECTRAL_UPPER - SPECTRAL_LOWER) / 2
H = build_hamiltonian()
scaled_H = H / halfwidth
coefficients = projector_coefficients(MOMENTS, (exact_mu - center) / halfwidth)

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
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "diagnostic" => "coordinate_interleaved_hadamard_vs_saved_ed",
        "reference_directory" => reference_directory,
        "backend" => string(backend),
        "device_name" => device_name,
        "device_total_memory_bytes" => device_total_memory,
        "device_free_memory_before_bytes" => device_free_memory_before,
        "nx" => NX,
        "ny" => NY,
        "matrix_dimension" => N,
        "moments" => MOMENTS,
        "spectral_lower" => SPECTRAL_LOWER,
        "spectral_upper" => SPECTRAL_UPPER,
        "probe_ordering" => "interleaved_coordinate_bits",
        "created_at" => string(now(UTC)),
    ))
end

summary_path = joinpath(output, "summary.csv")
open(summary_path, "w") do io
    write_csv_row(io, (
        "probes", "seed", "kpm_time_s", "measurement_time_s",
        "trace_error", "density_max_abs_error", "density_rms_error",
        "horizontal_bond_max_abs_error", "horizontal_bond_rms_error",
        "vertical_bond_max_abs_error", "vertical_bond_rms_error",
    ))
end

for (probes, seed) in PROBE_CONFIGURATIONS
    println("Coordinate-Hadamard KPM: R=$probes seed=$seed")
    flush(stdout)
    host_probes = coded_hadamard_matrix(codes, probes, seed)
    probes_backend = backend == :cuda ? CUDA.CuArray(host_probes) : host_probes
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
    filtered = backend == :cuda ? Array(filtered_backend) : filtered_backend
    measurement = @timed local_estimates(filtered, host_probes, bonds)
    density = measurement.value.density
    bond_order = measurement.value.bond_order
    density_error = density - exact_density
    bond_error = bond_order - exact_bond_order

    open(summary_path, "a") do io
        write_csv_row(io, (
            probes,
            seed,
            kpm_calculation.time,
            measurement.time,
            abs(sum(density) - NE),
            maximum(abs, density_error),
            norm(density_error) / sqrt(N),
            maximum(abs, @view bond_error[horizontal]),
            norm(@view bond_error[horizontal]) / sqrt(length(horizontal)),
            maximum(abs, @view bond_error[vertical]),
            norm(@view bond_error[vertical]) / sqrt(length(vertical)),
        ))
    end
    @printf("  density max/RMS error = %.6e / %.6e\n",
        maximum(abs, density_error), norm(density_error) / sqrt(N))
    flush(stdout)

    host_probes = nothing
    probes_backend = nothing
    filtered_backend = nothing
    filtered = nothing
    GC.gc()
    backend == :cuda && CUDA.reclaim()
end

println("Coordinate-Hadamard comparison complete: $output")

