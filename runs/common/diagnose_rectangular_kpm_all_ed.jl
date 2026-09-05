#!/usr/bin/env julia

"""All-site/all-bond KPM validation against one reusable dense ED reference.

The fixed open rectangle has 128x64=8192 sites. Dense diagonalization is
performed once. Several Hadamard probe counts are then evaluated on the GPU,
and every site density and every horizontal/vertical nearest-neighbour
density-matrix element is compared with the occupied eigenspace.
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 1:2 || error(
    "usage: diagnose_rectangular_kpm_all_ed.jl OUTPUT_DIRECTORY [BACKEND]",
)
output = abspath(ARGS[1])
backend = length(ARGS) == 2 ? Symbol(lowercase(ARGS[2])) : :cpu
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
    (512, 20260730),
    (1024, 20260730),
    (2048, 20260730),
]

site_index(x, y) = x + NX * y + 1
tx(x) = -1.0 - V2 * cos(2π * TAU * (Float64(x) + 0.5))
ty(y) = -1.0 - V2 * cos(2π * TAU * (Float64(y) + 0.5))
seed_field(x, y) = iseven(x + y) ? SEED_AMPLITUDE : -SEED_AMPLITUDE

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))

function build_hamiltonian_and_bonds()
    rows = Int[]
    columns = Int[]
    values = Float64[]
    bonds = Tuple{Int,Int,Symbol}[]
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
            push!(bonds, (site, neighbour, :horizontal))
        end
        if y < NY - 1
            neighbour = site_index(x, y + 1)
            hopping = ty(y)
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

function one_body_energy(H, density, bonds, bond_order)
    onsite = sum(H[site, site] * density[site] for site in eachindex(density))
    hopping = sum(
        2H[site, neighbour] * bond_order[index]
        for (index, (site, neighbour, _)) in enumerate(bonds)
    )
    return onsite + hopping
end

function checkerboard_order(density)
    return sum(0:(NY - 1)) do y
        sum(0:(NX - 1)) do x
            (iseven(x + y) ? 1.0 : -1.0) * density[site_index(x, y)]
        end
    end / N
end

mkpath(output)
H, bonds = build_hamiltonian_and_bonds()
horizontal = findall(bond -> bond[3] == :horizontal, bonds)
vertical = findall(bond -> bond[3] == :vertical, bonds)

println("Dense reference: rectangle=$(NX)x$(NY), N=$N, Ne=$NE")
flush(stdout)
exact_calculation = @timed begin
    eigenpairs = eigen(Symmetric(Matrix(H)))
    occupied = @view eigenpairs.vectors[:, 1:NE]
    density = vec(sum(abs2, occupied; dims=2))
    bond_order = [
        dot(@view(occupied[site, :]), @view(occupied[neighbour, :]))
        for (site, neighbour, _) in bonds
    ]
    (; eigenpairs, density, bond_order)
end
eigenpairs = exact_calculation.value.eigenpairs
exact_density = exact_calculation.value.density
exact_bond_order = exact_calculation.value.bond_order
exact_eigenvalues = eigenpairs.values
lambda_min = first(exact_eigenvalues)
lambda_max = last(exact_eigenvalues)
SPECTRAL_LOWER <= lambda_min && SPECTRAL_UPPER >= lambda_max || error(
    "safe interval [$SPECTRAL_LOWER,$SPECTRAL_UPPER] does not contain " *
    "the exact spectrum [$lambda_min,$lambda_max]",
)
exact_mu = (exact_eigenvalues[NE] + exact_eigenvalues[NE + 1]) / 2
exact_energy = sum(@view exact_eigenvalues[1:NE])
exact_diagonalization_time = exact_calculation.time
exact_diagonalization_bytes = exact_calculation.bytes
# Only local observables and eigenvalues are needed below. Release the dense
# eigenvector matrix before allocating the largest probe blocks.
eigenpairs = nothing
exact_calculation = nothing
GC.gc()

center = (SPECTRAL_LOWER + SPECTRAL_UPPER) / 2
halfwidth = (SPECTRAL_UPPER - SPECTRAL_LOWER) / 2
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, N, N)) / halfwidth
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

metadata = Dict(
    "diagnostic" => "rectangular_all_local_kpm_vs_ed",
    "backend" => string(backend),
    "device_name" => device_name,
    "device_total_memory_bytes" => device_total_memory,
    "device_free_memory_before_bytes" => device_free_memory_before,
    "nx" => NX,
    "ny" => NY,
    "matrix_dimension" => N,
    "target_particles" => NE,
    "V2" => V2,
    "tau" => TAU,
    "seed_amplitude" => SEED_AMPLITUDE,
    "moments" => MOMENTS,
    "spectral_lower" => SPECTRAL_LOWER,
    "spectral_upper" => SPECTRAL_UPPER,
    "exact_lambda_min" => lambda_min,
    "exact_lambda_max" => lambda_max,
    "exact_homo" => exact_eigenvalues[NE],
    "exact_lumo" => exact_eigenvalues[NE + 1],
    "exact_fermi_gap" => exact_eigenvalues[NE + 1] - exact_eigenvalues[NE],
    "exact_diagonalization_time_s" => exact_diagonalization_time,
    "exact_diagonalization_allocations_bytes" => exact_diagonalization_bytes,
    "probe_configurations" => [
        "R=$(configuration[1]),seed=$(configuration[2])"
        for configuration in PROBE_CONFIGURATIONS
    ],
    "created_at" => string(now(UTC)),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

summary_path = joinpath(output, "summary.csv")
open(summary_path, "w") do io
    write_csv_row(io, (
        "probes", "seed", "moments", "kpm_time_s",
        "transfer_to_device_time_s", "transfer_to_host_time_s",
        "measurement_time_s", "trace", "trace_error",
        "density_max_abs_error", "density_rms_error",
        "horizontal_bond_max_abs_error", "horizontal_bond_rms_error",
        "vertical_bond_max_abs_error", "vertical_bond_rms_error",
        "checkerboard_order_error", "one_body_energy_error",
    ))
end

for (probes, seed) in PROBE_CONFIGURATIONS
    case_name = "m$(MOMENTS)_hadamard_r$(probes)_seed$(seed)"
    case_output = joinpath(output, case_name)
    mkpath(case_output)
    println("KPM case: $case_name")
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
        filtered_backend = kpm_apply(
            H_backend,
            probes_backend,
            coefficients;
            synchronize=synchronize,
        )
        synchronize()
        filtered_backend
    end
    filtered_backend = kpm_calculation.value
    transfer_to_host = @timed begin
        filtered = backend == :cuda ? Array(filtered_backend) :
            filtered_backend
        filtered
    end
    filtered = transfer_to_host.value

    measurement = @timed local_estimates(filtered, host_probes, bonds)
    density = measurement.value.density
    bond_order = measurement.value.bond_order
    density_error = density - exact_density
    bond_error = bond_order - exact_bond_order
    energy = one_body_energy(H, density, bonds, bond_order)

    values = (
        probes,
        seed,
        MOMENTS,
        kpm_calculation.time,
        transfer_to_device.time,
        transfer_to_host.time,
        measurement.time,
        sum(density),
        abs(sum(density) - NE),
        maximum(abs, density_error),
        norm(density_error) / sqrt(N),
        maximum(abs, @view bond_error[horizontal]),
        norm(@view bond_error[horizontal]) / sqrt(length(horizontal)),
        maximum(abs, @view bond_error[vertical]),
        norm(@view bond_error[vertical]) / sqrt(length(vertical)),
        checkerboard_order(density) - checkerboard_order(exact_density),
        energy - exact_energy,
    )
    open(summary_path, "a") do io
        write_csv_row(io, values)
    end

    open(joinpath(case_output, "density.csv"), "w") do io
        write_csv_row(io, ("site", "x", "y", "exact", "kpm", "error"))
        for y in 0:(NY - 1), x in 0:(NX - 1)
            site = site_index(x, y)
            write_csv_row(io, (
                site, x, y, exact_density[site], density[site],
                density_error[site],
            ))
        end
    end
    open(joinpath(case_output, "bonds.csv"), "w") do io
        write_csv_row(io, (
            "site", "neighbour", "orientation", "exact", "kpm", "error",
        ))
        for (index, bond) in enumerate(bonds)
            write_csv_row(io, (
                bond[1], bond[2], bond[3], exact_bond_order[index],
                bond_order[index], bond_error[index],
            ))
        end
    end

    @printf("  density max/RMS error = %.6e / %.6e\n",
        maximum(abs, density_error), norm(density_error) / sqrt(N))
    flush(stdout)

    # Release the largest per-case arrays before constructing the next probe
    # block. CUDA's pool may retain reserved memory for efficient reuse.
    host_probes = nothing
    probes_backend = nothing
    filtered_backend = nothing
    filtered = nothing
    GC.gc()
    backend == :cuda && CUDA.reclaim()
end

println("All-site rectangular KPM diagnostic complete: $output")
