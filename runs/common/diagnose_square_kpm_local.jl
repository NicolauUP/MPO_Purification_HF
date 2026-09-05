#!/usr/bin/env julia

"""Validate sparse-vector KPM local observables against exact diagonalization.

This diagnostic deliberately does not construct an MPO projector and does not
run SCF. It approximates the fixed-Hamiltonian zero-temperature projector

    rho = Theta(mu I - H)

with a Jackson-damped Chebyshev expansion, applies that polynomial to a block
of probing vectors, and estimates only the density and open-boundary
nearest-neighbour density-matrix elements.

A complete Hadamard basis (`PROBE_METHOD=hadamard`, `PROBES=N`) removes probing
error exactly and therefore isolates the polynomial approximation error.
Smaller randomized Hadamard or Rademacher sets include probing error.
"""

using Dates
using LinearAlgebra
using MPO_MeanField
using Printf
using Random
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) in 8:9 || error(
    "usage: diagnose_square_kpm_local.jl CAMPAIGN_FILE TASK_INDEX MOMENTS " *
    "PROBE_METHOD PROBES SEED OUTPUT_DIRECTORY BACKEND [SPECTRAL_BOUND]",
)

campaign_file = abspath(ARGS[1])
task_index = parse(Int, ARGS[2])
moments = parse(Int, ARGS[3])
probe_method = Symbol(lowercase(ARGS[4]))
probes = parse(Int, ARGS[5])
seed = parse(Int, ARGS[6])
output = abspath(ARGS[7])
backend = Symbol(lowercase(ARGS[8]))
spectral_bound = length(ARGS) == 9 ? parse(Float64, ARGS[9]) : nothing

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
moments >= 2 || error("MOMENTS must be at least 2")
probe_method in (:hadamard, :rademacher) ||
    error("PROBE_METHOD must be hadamard or rademacher")
probes > 0 || error("PROBES must be positive")
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
            append!(columns, (neighbour,site))
            append!(values, (hopping, hopping))
        end
        if y < side - 1
            neighbour = square_lattice_index(x, y + 1, params.L)
            hopping = ty(x, y)
            append!(rows, (site, neighbour))
            append!(columns, (neighbour,site))
            append!(values, (hopping, hopping))
        end
    end
    return sparse(rows, columns, values, N, N)
end

function local_estimates(filtered, probes, bonds)
    probe_count = size(probes, 2)
    density_samples = filtered .* probes
    density = vec(mean(density_samples; dims=2))
    density_variance = probe_count > 1 ?
        vec(var(density_samples; dims=2, corrected=true)) :
        zeros(eltype(density), length(density))
    bond_order = Vector{Float64}(undef, length(bonds))
    bond_samples = Matrix{Float64}(undef, length(bonds), probe_count)
    for (index, (site, neighbour, _)) in enumerate(bonds)
        # Symmetrization reduces finite-probe noise and enforces the real
        # Hermitian convention used by this Hamiltonian.
        @views bond_samples[index, :] .= (
            filtered[site, :] .* probes[neighbour, :] .+
            filtered[neighbour, :] .* probes[site, :]
        ) ./ 2
        bond_order[index] = mean(@view bond_samples[index, :])
    end
    bond_variance = probe_count > 1 ?
        vec(var(bond_samples; dims=2, corrected=true)) :
        zeros(eltype(bond_order), length(bond_order))
    return (;
        density,
        density_variance,
        density_samples,
        bond_order,
        bond_variance,
        bond_samples,
    )
end

relative_standard_error(variance, average, probe_count) =
    sqrt(variance / probe_count) / abs(average)
coefficient_of_variation(variance, average) =
    sqrt(variance) / abs(average)
variance_to_abs_mean(variance, average) = variance / abs(average)

function global_relative_standard_error(variance, estimate, probe_count)
    return sqrt(mean(variance) / probe_count) / sqrt(mean(abs2, estimate))
end

function checkerboard_order(density, L)
    return sum(eachindex(density)) do site
        x, y = square_lattice_decoder(site - 1, L)
        (iseven(x + y) ? 1.0 : -1.0) * density[site]
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

N = 2^params.L
probes <= N || error("PROBES=$probes exceeds N=$N")
Ne = round(Int, params.density * N)
bonds = collect(square_undirected_bonds(params.L))
started_at = now(UTC)
mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))

H = direct_initial_hamiltonian(params)
exact_calculation = @timed begin
    eigenpairs = eigen(Symmetric(Matrix(H)))
    occupied = @view eigenpairs.vectors[:, 1:Ne]
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
exact_mu = (eigenpairs.values[Ne] + eigenpairs.values[Ne + 1]) / 2

lower, upper = if isnothing(spectral_bound)
    validate_spectral_bounds(spec.spectral_bounds...)
else
    (-spectral_bound, spectral_bound)
end
lambda_min = first(eigenpairs.values)
lambda_max = last(eigenpairs.values)
lower <= lambda_min && upper >= lambda_max || error(
    "supplied KPM interval [$lower,$upper] does not contain " *
    "the exact spectrum [$lambda_min,$lambda_max]",
)
center = (lower + upper) / 2
halfwidth = (upper - lower) / 2
scaled_mu = (exact_mu - center) / halfwidth
scaled_H = (H - center * sparse(I, N, N)) / halfwidth
coefficients = projector_coefficients(moments, scaled_mu)
host_probes = probing_matrix(N, probes, probe_method, seed)

device_name = "CPU"
device_total_memory = 0
device_free_memory_before = 0
transfer_to_device_time = 0.0
synchronize = () -> nothing
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    device_total_memory = CUDA.total_memory()
    device_free_memory_before = CUDA.free_memory()
    transfer = @timed begin
        scaled_H_backend = CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H)
        probes_backend = CUDA.CuArray(host_probes)
        (; scaled_H_backend, probes_backend)
    end
    transfer_to_device_time = transfer.time
    scaled_H_backend = transfer.value.scaled_H_backend
    probes_backend = transfer.value.probes_backend
    synchronize = CUDA.synchronize
else
    scaled_H_backend = scaled_H
    probes_backend = host_probes
end

synchronize()
kpm_calculation = @timed begin
    filtered_backend = kpm_apply(
        scaled_H_backend,
        probes_backend,
        coefficients;
        synchronize=synchronize,
    )
    synchronize()
    filtered_backend
end
filtered = backend == :cuda ? Array(kpm_calculation.value) :
    kpm_calculation.value

measurement = @timed local_estimates(filtered, host_probes, bonds)
kpm_density = measurement.value.density
density_variance = measurement.value.density_variance
density_samples = measurement.value.density_samples
kpm_bond_order = measurement.value.bond_order
bond_variance = measurement.value.bond_variance
bond_samples = measurement.value.bond_samples
density_error = kpm_density - exact_density
bond_error = kpm_bond_order - exact_bond_order
horizontal = findall(bond -> bond[3] == :horizontal, bonds)
vertical = findall(bond -> bond[3] == :vertical, bonds)
exact_energy = sum(@view eigenpairs.values[1:Ne])
kpm_energy = one_body_energy(H, kpm_density, bonds, kpm_bond_order)

summary = Dict(
    "campaign" => campaign.name,
    "label" => spec.label,
    "task_index" => task_index,
    "diagnostic" => "sparse_vector_jackson_kpm_local",
    "backend" => string(backend),
    "device_name" => device_name,
    "device_total_memory_bytes" => device_total_memory,
    "device_free_memory_before_bytes" => device_free_memory_before,
    "matrix_dimension" => N,
    "target_particles" => Ne,
    "moments" => moments,
    "probe_method" => string(probe_method),
    "probes" => probes,
    "seed" => seed,
    "spectral_lower" => lower,
    "spectral_upper" => upper,
    "exact_lambda_min" => lambda_min,
    "exact_lambda_max" => lambda_max,
    "exact_homo" => eigenpairs.values[Ne],
    "exact_lumo" => eigenpairs.values[Ne + 1],
    "exact_fermi_gap" => eigenpairs.values[Ne + 1] - eigenpairs.values[Ne],
    "exact_chemical_potential" => exact_mu,
    "trace" => sum(kpm_density),
    "trace_error" => abs(sum(kpm_density) - Ne),
    "density_max_abs_error" => maximum(abs, density_error),
    "density_rms_error" => norm(density_error) / sqrt(N),
    "density_probe_naive_relative_standard_error" =>
        global_relative_standard_error(density_variance, kpm_density, probes),
    "probing_reconstruction_exact" =>
        probe_method == :hadamard && probes == N,
    "horizontal_bond_max_abs_error" =>
        maximum(abs, @view bond_error[horizontal]),
    "horizontal_bond_rms_error" =>
        norm(@view bond_error[horizontal]) / sqrt(length(horizontal)),
    "horizontal_bond_probe_naive_relative_standard_error" =>
        global_relative_standard_error(
            @view(bond_variance[horizontal]),
            @view(kpm_bond_order[horizontal]),
            probes,
        ),
    "vertical_bond_max_abs_error" =>
        maximum(abs, @view bond_error[vertical]),
    "vertical_bond_rms_error" =>
        norm(@view bond_error[vertical]) / sqrt(length(vertical)),
    "vertical_bond_probe_naive_relative_standard_error" =>
        global_relative_standard_error(
            @view(bond_variance[vertical]),
            @view(kpm_bond_order[vertical]),
            probes,
        ),
    "exact_checkerboard_order" => checkerboard_order(exact_density, params.L),
    "kpm_checkerboard_order" => checkerboard_order(kpm_density, params.L),
    "exact_one_body_energy" => exact_energy,
    "kpm_one_body_energy" => kpm_energy,
    "one_body_energy_error" => kpm_energy - exact_energy,
    "exact_diagonalization_time_s" => exact_calculation.time,
    "exact_diagonalization_allocations_bytes" => exact_calculation.bytes,
    "transfer_to_device_time_s" => transfer_to_device_time,
    "kpm_time_s" => kpm_calculation.time,
    "kpm_allocations_bytes" => kpm_calculation.bytes,
    "measurement_time_s" => measurement.time,
    "measurement_allocations_bytes" => measurement.bytes,
    "started_at" => string(started_at),
    "finished_at" => string(now(UTC)),
)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, summary)
end

open(joinpath(output, "density.csv"), "w") do io
    write_csv_row(io, (
        "site", "x", "y", "exact", "kpm", "error",
        "probe_variance", "variance_to_abs_mean",
        "coefficient_of_variation", "naive_relative_standard_error",
    ))
    for site in 1:N
        x, y = square_lattice_decoder(site - 1, params.L)
        average = kpm_density[site]
        variance = density_variance[site]
        write_csv_row(io, (
            site, x, y, exact_density[site], average, density_error[site],
            variance, variance_to_abs_mean(variance, average),
            coefficient_of_variation(variance, average),
            relative_standard_error(variance, average, probes),
        ))
    end
end
open(joinpath(output, "bonds.csv"), "w") do io
    write_csv_row(io, (
        "site", "neighbour", "orientation", "exact", "kpm", "error",
        "probe_variance", "variance_to_abs_mean",
        "coefficient_of_variation", "naive_relative_standard_error",
    ))
    for (index, bond) in enumerate(bonds)
        average = kpm_bond_order[index]
        variance = bond_variance[index]
        write_csv_row(io, (
            bond[1], bond[2], bond[3], exact_bond_order[index], average,
            bond_error[index], variance,
            variance_to_abs_mean(variance, average),
            coefficient_of_variation(variance, average),
            relative_standard_error(variance, average, probes),
        ))
    end
end

open(joinpath(output, "probe_statistics.csv"), "w") do io
    write_csv_row(io, (
        "probe", "density_sample_average", "density_error_variance",
        "density_error_variance_to_exact_average",
        "density_relative_rms_error", "bond_sample_average",
        "bond_error_variance", "bond_error_variance_to_exact_abs_average",
        "bond_relative_rms_error",
    ))
    exact_density_average = mean(exact_density)
    exact_bond_abs_average = mean(abs, exact_bond_order)
    exact_density_rms = sqrt(mean(abs2, exact_density))
    exact_bond_rms = sqrt(mean(abs2, exact_bond_order))
    for probe in 1:probes
        density_probe_error =
            @view(density_samples[:, probe]) .- exact_density
        bond_probe_error =
            @view(bond_samples[:, probe]) .- exact_bond_order
        write_csv_row(io, (
            probe,
            mean(@view density_samples[:, probe]),
            var(density_probe_error; corrected=false),
            var(density_probe_error; corrected=false) /
                abs(exact_density_average),
            sqrt(mean(abs2, density_probe_error)) / exact_density_rms,
            mean(@view bond_samples[:, probe]),
            var(bond_probe_error; corrected=false),
            var(bond_probe_error; corrected=false) /
                exact_bond_abs_average,
            sqrt(mean(abs2, bond_probe_error)) / exact_bond_rms,
        ))
    end
end

println("KPM diagnostic complete: $output")
println("  N=$N moments=$moments probes=$probes method=$probe_method")
@printf("  density max/RMS error = %.6e / %.6e\n",
    summary["density_max_abs_error"], summary["density_rms_error"])
@printf("  trace error = %.6e, energy error = %.6e\n",
    summary["trace_error"], summary["one_body_energy_error"])
