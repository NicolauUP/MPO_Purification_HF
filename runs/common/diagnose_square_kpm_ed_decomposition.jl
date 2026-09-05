#!/usr/bin/env julia

"""Decompose KPM polynomial and probing errors against exact diagonalization.

For a fixed square Hamiltonian, ED supplies both the exact occupied projector
and deterministic local observables of each finite Jackson--Chebyshev
polynomial. KPM probing is then compared with both references:

    total error      = probed KPM - exact projector
    probing error    = probed KPM - deterministic KPM polynomial
    polynomial error = deterministic KPM polynomial - exact projector

No dense projector matrix is formed.
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
    "usage: diagnose_square_kpm_ed_decomposition.jl " *
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

const MOMENT_COUNTS = (400, 800, 1200, 1600)
const PRIMARY_SEED = 510578
const PROBE_CONFIGURATIONS = (
    (512, PRIMARY_SEED),
    (1024, PRIMARY_SEED),
    (2048, PRIMARY_SEED),
    (4096, PRIMARY_SEED),
    (1024, 20260730),
    (2048, 20260730),
    (4096, 20260730),
)

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

function polynomial_values(coefficients, scaled_eigenvalues)
    degree = length(coefficients) - 1
    previous = ones(Float64, length(scaled_eigenvalues))
    current = copy(scaled_eigenvalues)
    result = (coefficients[1] / 2) .* previous .+
             coefficients[2] .* current
    for order in 2:degree
        following = 2 .* scaled_eigenvalues .* current .- previous
        result .+= coefficients[order + 1] .* following
        previous, current = current, following
    end
    return result
end

function spectral_local_observables(site_spectra, weights, bonds)
    # site_spectra[:, site] contains all eigenstate amplitudes at one site.
    # This transpose of the usual eigenvector layout makes the threaded local
    # contractions contiguous in memory.
    N = size(site_spectra, 2)
    density = Vector{Float64}(undef, N)
    Threads.@threads for site in 1:N
        value = 0.0
        @inbounds @simd for state in axes(site_spectra, 1)
            amplitude = site_spectra[state, site]
            value += weights[state] * amplitude * amplitude
        end
        density[site] = value
    end

    bond_order = Vector{Float64}(undef, length(bonds))
    Threads.@threads for index in eachindex(bonds)
        site, neighbour, _ = bonds[index]
        value = 0.0
        @inbounds @simd for state in axes(site_spectra, 1)
            value += weights[state] * site_spectra[state, site] *
                     site_spectra[state, neighbour]
        end
        bond_order[index] = value
    end
    return (; density, bond_order)
end

function probed_local_observables(filtered, probes, bonds)
    R = size(probes, 2)
    density = vec(mean(filtered .* probes; dims=2))
    bond_order = Vector{Float64}(undef, length(bonds))
    Threads.@threads for index in eachindex(bonds)
        site, neighbour, _ = bonds[index]
        bond_order[index] = (
            dot(@view(filtered[site, :]), @view(probes[neighbour, :])) +
            dot(@view(filtered[neighbour, :]), @view(probes[site, :]))
        ) / (2R)
    end
    return (; density, bond_order)
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

function errors(values, reference)
    difference = values - reference
    return (
        maximum(abs, difference),
        norm(difference) / sqrt(length(difference)),
    )
end

N = 2^params.L
side = 2^div(params.L, 2)
N == 16_384 || error("this validation expects N=16,384, got N=$N")
Ne = round(Int, params.density * N)
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
center = (lower + upper) / 2
halfwidth = (upper - lower) / 2
H, bonds = direct_initial_hamiltonian_and_bonds(params)
horizontal = findall(bond -> bond[3] == :horizontal, bonds)
vertical = findall(bond -> bond[3] == :vertical, bonds)

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
println("Exact diagonalization: side=$side N=$N Ne=$Ne")
flush(stdout)
diagonalization = @timed eigen!(Symmetric(Matrix(H)))
diagonalization_time = diagonalization.time
diagonalization_bytes = diagonalization.bytes
eigenpairs = diagonalization.value
eigenvalues = copy(eigenpairs.values)
# Eigenvectors are normally indexed [site, state]. Transpose once so that all
# states contributing to one local observable are contiguous.
site_spectra = Matrix(transpose(eigenpairs.vectors))
diagonalization = nothing
eigenpairs = nothing
GC.gc()
lower <= first(eigenvalues) && last(eigenvalues) <= upper || error(
    "exact spectrum [$(first(eigenvalues)), $(last(eigenvalues))] " *
    "lies outside configured bounds [$lower, $upper]",
)
homo = eigenvalues[Ne]
lumo = eigenvalues[Ne + 1]
mu = (homo + lumo) / 2
scaled_eigenvalues = (eigenvalues .- center) ./ halfwidth

println("Exact local projector observables")
flush(stdout)
exact_measurement = @timed spectral_local_observables(
    site_spectra,
    vcat(ones(Float64, Ne), zeros(Float64, N - Ne)),
    bonds,
)
exact = exact_measurement.value
exact_order = checkerboard_order(exact.density, params)
exact_energy = one_body_energy(H, exact.density, bonds, exact.bond_order)

polynomial_references = Dict{Int,NamedTuple}()
polynomial_path = joinpath(output, "polynomial_errors.csv")
open(polynomial_path, "w") do io
    write_csv_row(io, (
        "moments", "trace", "trace_error", "occupation_min",
        "occupation_max", "density_max_abs_error", "density_rms_error",
        "horizontal_bond_max_abs_error", "horizontal_bond_rms_error",
        "vertical_bond_max_abs_error", "vertical_bond_rms_error",
        "checkerboard_order_error", "one_body_energy_error",
        "local_measurement_time_s",
    ))
end

for moments in MOMENT_COUNTS
    println("Deterministic polynomial reference: M=$moments")
    flush(stdout)
    coefficients = projector_coefficients(moments, (mu - center) / halfwidth)
    occupations = polynomial_values(coefficients, scaled_eigenvalues)
    measurement = @timed spectral_local_observables(
        site_spectra, occupations, bonds,
    )
    value = measurement.value
    polynomial_references[moments] = (;
        coefficients, density=value.density, bond_order=value.bond_order,
    )
    density_error = errors(value.density, exact.density)
    horizontal_error = errors(
        @view(value.bond_order[horizontal]),
        @view(exact.bond_order[horizontal]),
    )
    vertical_error = errors(
        @view(value.bond_order[vertical]),
        @view(exact.bond_order[vertical]),
    )
    open(polynomial_path, "a") do io
        write_csv_row(io, (
            moments, sum(value.density), abs(sum(value.density) - Ne),
            minimum(occupations), maximum(occupations), density_error...,
            horizontal_error..., vertical_error...,
            checkerboard_order(value.density, params) - exact_order,
            one_body_energy(H, value.density, bonds, value.bond_order) -
                exact_energy,
            measurement.time,
        ))
    end
end

device_name = "CPU"
device_total_memory = 0
device_free_memory_before = 0
synchronize = () -> nothing
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, N, N)) / halfwidth
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    device_total_memory = CUDA.total_memory()
    device_free_memory_before = CUDA.free_memory()
    H_backend = CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H)
    synchronize = CUDA.synchronize
else
    H_backend = scaled_H
end

summary_path = joinpath(output, "probing_errors.csv")
open(summary_path, "w") do io
    write_csv_row(io, (
        "moments", "probes", "seed", "kpm_time_s",
        "transfer_to_device_time_s", "transfer_to_host_time_s",
        "measurement_time_s", "trace", "trace_error",
        "total_density_max_abs_error", "total_density_rms_error",
        "probing_density_max_abs_error", "probing_density_rms_error",
        "total_horizontal_bond_max_abs_error",
        "total_horizontal_bond_rms_error",
        "probing_horizontal_bond_max_abs_error",
        "probing_horizontal_bond_rms_error",
        "total_vertical_bond_max_abs_error",
        "total_vertical_bond_rms_error",
        "probing_vertical_bond_max_abs_error",
        "probing_vertical_bond_rms_error",
        "checkerboard_order_error", "one_body_energy_error",
    ))
end

for moments in MOMENT_COUNTS
    polynomial = polynomial_references[moments]
    for (probes, seed) in PROBE_CONFIGURATIONS
        println("Probed KPM: M=$moments R=$probes seed=$seed")
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
                H_backend, probes_backend, polynomial.coefficients;
                synchronize=synchronize,
            )
            synchronize()
            value
        end
        filtered_backend = kpm_calculation.value
        transfer_to_host = @timed begin
            value = backend == :cuda ? Array(filtered_backend) :
                filtered_backend
            value
        end
        filtered = transfer_to_host.value
        measurement = @timed probed_local_observables(
            filtered, host_probes, bonds,
        )
        value = measurement.value

        total_density = errors(value.density, exact.density)
        probing_density = errors(value.density, polynomial.density)
        total_horizontal = errors(
            @view(value.bond_order[horizontal]),
            @view(exact.bond_order[horizontal]),
        )
        probing_horizontal = errors(
            @view(value.bond_order[horizontal]),
            @view(polynomial.bond_order[horizontal]),
        )
        total_vertical = errors(
            @view(value.bond_order[vertical]),
            @view(exact.bond_order[vertical]),
        )
        probing_vertical = errors(
            @view(value.bond_order[vertical]),
            @view(polynomial.bond_order[vertical]),
        )
        open(summary_path, "a") do io
            write_csv_row(io, (
                moments, probes, seed, kpm_calculation.time,
                transfer_to_device.time, transfer_to_host.time,
                measurement.time, sum(value.density),
                abs(sum(value.density) - Ne), total_density...,
                probing_density..., total_horizontal..., probing_horizontal...,
                total_vertical..., probing_vertical...,
                checkerboard_order(value.density, params) - exact_order,
                one_body_energy(H, value.density, bonds, value.bond_order) -
                    exact_energy,
            ))
        end

        host_probes = nothing
        probes_backend = nothing
        filtered_backend = nothing
        filtered = nothing
        GC.gc()
        backend == :cuda && CUDA.reclaim()
    end
end

metadata = Dict(
    "campaign" => campaign.name,
    "label" => spec.label,
    "diagnostic" => "square_kpm_ed_error_decomposition",
    "backend" => string(backend),
    "device_name" => device_name,
    "device_total_memory_bytes" => device_total_memory,
    "device_free_memory_before_bytes" => device_free_memory_before,
    "side" => side,
    "matrix_dimension" => N,
    "target_particles" => Ne,
    "spectral_lower" => lower,
    "spectral_upper" => upper,
    "exact_lambda_min" => first(eigenvalues),
    "exact_lambda_max" => last(eigenvalues),
    "exact_homo" => homo,
    "exact_lumo" => lumo,
    "exact_fermi_gap" => lumo - homo,
    "chemical_potential" => mu,
    "exact_checkerboard_order" => exact_order,
    "exact_one_body_energy" => exact_energy,
    "moment_counts" => collect(MOMENT_COUNTS),
    "probe_configurations" => [
        "R=$(configuration[1]),seed=$(configuration[2])"
        for configuration in PROBE_CONFIGURATIONS
    ],
    "diagonalization_time_s" => diagonalization_time,
    "diagonalization_allocations_bytes" => diagonalization_bytes,
    "exact_local_measurement_time_s" => exact_measurement.time,
    "julia_threads" => Threads.nthreads(),
    "finished_at" => string(now(UTC)),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

println("KPM/ED decomposition complete: $output")
