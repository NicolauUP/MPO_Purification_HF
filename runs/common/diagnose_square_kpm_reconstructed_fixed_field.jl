#!/usr/bin/env julia

"""Separate KPM polynomial and probing errors on a saved square-HF field.

The production square runner stores local observables but not the mixed
Hartree/Fock vectors themselves. This diagnostic therefore reconstructs a
fixed effective Hamiltonian from the *audited* saved density and bond order,
then holds that reconstructed field fixed. It does not run SCF. Consequently,
it diagnoses the KPM representation of a near-self-consistent field; it does
not replace a final SCF-stationarity check.
"""

using Dates
using LinearAlgebra
using MPO_MeanField
using Printf
using SparseArrays
using Statistics
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) == 5 || error(
    "usage: diagnose_square_kpm_reconstructed_fixed_field.jl " *
    "CAMPAIGN_FILE CASE_INDEX SOURCE_RESULT OUTPUT_DIRECTORY BACKEND",
)
campaign_file = abspath(ARGS[1])
case_index = parse(Int, ARGS[2])
source_result = abspath(ARGS[3])
output = abspath(ARGS[4])
backend = Symbol(lowercase(ARGS[5]))

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
isdir(source_result) || error("source result does not exist: $source_result")
backend in (:cpu, :cuda) || error("BACKEND must be cpu or cuda")
ispath(output) && error("refusing to overwrite existing output directory: $output")

# The seed only selects a campaign run label; it does not enter the
# reconstructed Hamiltonian, whose onsite term is W and whose fields are read
# from SOURCE_RESULT.
ENV["KPM_SEED_STABILITY_CASE"] = string(case_index)
ENV["KPM_SEED_STABILITY_SEED"] = "1"
include(campaign_file)
@isdefined(campaign) || error("campaign file must define `campaign`")
spec = only(campaign.runs)
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")

if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
end

const CONFIGURATIONS = (
    (label="M3200_R2048", moments=3200, probes=2048, seed=20260730),
    (label="M3200_R4096", moments=3200, probes=4096, seed=20260730),
    (label="M4000_R2048", moments=4000, probes=2048, seed=20260730),
    (label="M4000_R4096", moments=4000, probes=4096, seed=20260730),
)

csv_escape(value) = '"' * replace(string(value), '"' => "\"\"") * '"'
write_csv_row(io, values) = println(io, join(csv_escape.(values), ','))
parse_csv_fields(line) = split(replace(chomp(line), '"' => ""), ',')

function lattice_data(params::ParametersSquare)
    side = 2^div(params.L, 2)
    N = side^2
    onsite = zeros(Float64, N)
    bonds = Tuple{Int,Int,Symbol}[]
    hopping = Float64[]
    tx(x, y) = params.t[1] isa Number ? Float64(params.t[1]) :
        Float64(params.t[1](x, y))
    ty(x, y) = params.t[2] isa Number ? Float64(params.t[2]) :
        Float64(params.t[2](x, y))
    for y in 0:(side - 1), x in 0:(side - 1)
        site = square_lattice_index(x, y, params.L)
        onsite[site] = isnothing(params.W) ? 0.0 : Float64(params.W(x, y))
        if x < side - 1
            push!(bonds, (site, square_lattice_index(x + 1, y, params.L), :horizontal))
            push!(hopping, tx(x, y))
        end
        if y < side - 1
            push!(bonds, (site, square_lattice_index(x, y + 1, params.L), :vertical))
            push!(hopping, ty(x, y))
        end
    end
    return (; side, N, onsite, bonds, hopping)
end

function read_site_column(path, expected_header, column, count)
    isfile(path) || error("missing $path")
    values = Vector{Float64}(undef, count)
    seen = falses(count)
    open(path) do io
        header = parse_csv_fields(readline(io))
        header == expected_header || error("unexpected header in $path: $header")
        site_index = only(findall(==("site"), header))
        column_index = only(findall(==(column), header))
        for line in eachline(io)
            row = parse_csv_fields(line)
            site = parse(Int, row[site_index])
            1 <= site <= count || error("invalid site $site in $path")
            seen[site] && error("duplicate site $site in $path")
            values[site] = parse(Float64, row[column_index])
            seen[site] = true
        end
    end
    all(seen) || error("missing site rows in $path")
    return values
end

function read_bond_column(path, expected_header, column, bonds)
    isfile(path) || error("missing $path")
    values = Vector{Float64}(undef, length(bonds))
    position = 0
    open(path) do io
        header = parse_csv_fields(readline(io))
        header == expected_header || error("unexpected header in $path: $header")
        site_index = only(findall(==("site"), header))
        neighbour_index = only(findall(==("neighbour"), header))
        orientation_index = only(findall(==("orientation"), header))
        column_index = only(findall(==(column), header))
        for line in eachline(io)
            position += 1
            position <= length(bonds) || error("too many bond rows in $path")
            row = parse_csv_fields(line)
            expected_site, expected_neighbour, expected_orientation = bonds[position]
            (parse(Int, row[site_index]), parse(Int, row[neighbour_index]), Symbol(row[orientation_index])) ==
                (expected_site, expected_neighbour, expected_orientation) ||
                error("bond ordering mismatch at row $position in $path")
            values[position] = parse(Float64, row[column_index])
        end
    end
    position == length(bonds) || error("expected $(length(bonds)) bond rows in $path, found $position")
    return values
end

function mean_fields(density, bond_order, data, U)
    hartree = zeros(Float64, data.N)
    for (site, neighbour, _) in data.bonds
        hartree[site] += U * density[neighbour]
        hartree[neighbour] += U * density[site]
    end
    return hartree, -U .* bond_order
end

function effective_hamiltonian(data, hartree, fock)
    rows = Vector{Int}(undef, data.N + 2length(data.bonds))
    columns = similar(rows)
    values = Vector{Float64}(undef, length(rows))
    position = 1
    for site in 1:data.N
        rows[position] = site; columns[position] = site
        values[position] = data.onsite[site] + hartree[site]
        position += 1
    end
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        coefficient = data.hopping[index] + fock[index]
        rows[position] = site; columns[position] = neighbour; values[position] = coefficient
        position += 1
        rows[position] = neighbour; columns[position] = site; values[position] = coefficient
        position += 1
    end
    return sparse(rows, columns, values, data.N, data.N)
end

function gershgorin_bounds(data, hartree, fock)
    radius = zeros(Float64, data.N)
    for index in eachindex(data.bonds)
        site, neighbour, _ = data.bonds[index]
        value = abs(data.hopping[index] + fock[index])
        radius[site] += value; radius[neighbour] += value
    end
    diagonal = data.onsite + hartree
    return minimum(diagonal - radius) - 0.1, maximum(diagonal + radius) + 0.1
end

function local_estimates(filtered, probes, bonds)
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

function checkerboard_order(density, data, L)
    return sum(0:(data.side - 1)) do y
        sum(0:(data.side - 1)) do x
            (iseven(x + y) ? 1.0 : -1.0) *
                density[square_lattice_index(x, y, L)]
        end
    end / data.N
end

function hf_energy(data, density, bond_order, U)
    kinetic = dot(data.onsite, density) + 2dot(data.hopping, bond_order)
    hartree = sum(data.bonds) do (site, neighbour, _)
        U * density[site] * density[neighbour]
    end
    fock = -U * sum(abs2, bond_order)
    return kinetic + hartree + fock
end

function metrics(values, reference)
    difference = values - reference
    return maximum(abs, difference), norm(difference) / sqrt(length(difference))
end

data = lattice_data(params)
N = data.N
Ne = round(Int, params.density * N)
density_path = joinpath(source_result, "density.csv")
bonds_path = joinpath(source_result, "bond_order.csv")
source_density = read_site_column(
    density_path, ["site", "x", "y", "production", "audit"], "audit", N,
)
source_bonds = read_bond_column(
    bonds_path, ["site", "neighbour", "orientation", "production", "audit"],
    "audit", data.bonds,
)
hartree, fock = mean_fields(source_density, source_bonds, data, Float64(params.U))
lower, upper = gershgorin_bounds(data, hartree, fock)
center, halfwidth = (lower + upper) / 2, (upper - lower) / 2
H = effective_hamiltonian(data, hartree, fock)
scaled_H = center == 0 ? H / halfwidth :
    (H - center * sparse(I, N, N)) / halfwidth

device_name = "CPU"
synchronize = () -> nothing
if backend == :cuda
    device_name = CUDA.name(CUDA.device())
    H_backend = CUDA.CUSPARSE.CuSparseMatrixCSR(scaled_H)
    synchronize = CUDA.synchronize
else
    H_backend = scaled_H
end

mkpath(output)
cp(campaign_file, joinpath(output, "input.jl"))
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "diagnostic" => "reconstructed_audited_fixed_field_kpm",
        "source_result" => source_result,
        "reconstruction" => "Hartree/Fock fields rebuilt from saved audit density and bond order",
        "backend" => string(backend),
        "device_name" => device_name,
        "matrix_dimension" => N,
        "target_particles" => Ne,
        "spectral_lower" => lower,
        "spectral_upper" => upper,
        "configurations" => [configuration.label for configuration in CONFIGURATIONS],
        "created_at" => string(now(UTC)),
    ))
end

results = Dict{String,NamedTuple}()
open(joinpath(output, "summary.csv"), "w") do io
    write_csv_row(io, (
        "configuration", "moments", "probes", "seed", "trace", "trace_error",
        "checkerboard_order", "energy_total", "kpm_time_s", "measurement_time_s",
        "source_density_max_abs_difference", "source_density_rms_difference",
        "source_bond_max_abs_difference", "source_bond_rms_difference",
    ))
    for configuration in CONFIGURATIONS
        println("Fixed-field KPM: $(configuration.label)")
        host_probes = probing_matrix(N, configuration.probes, :hadamard, configuration.seed)
        probes_backend = backend == :cuda ? CUDA.CuArray(host_probes) : host_probes
        synchronize()
        trace_moments = kpm_trace_moments(
            H_backend, probes_backend, configuration.moments; synchronize=synchronize,
        )
        mu_result = find_scaled_chemical_potential(trace_moments, Ne; tolerance=max(1e-6 * Ne, 1e-6))
        coefficients = projector_coefficients(configuration.moments, mu_result.scaled_mu)
        GC.gc(); backend == :cuda && CUDA.reclaim()
        kpm = @timed kpm_apply(
            H_backend, probes_backend, coefficients; synchronize=synchronize,
        )
        synchronize()
        filtered = backend == :cuda ? Array(kpm.value) : kpm.value
        measurement = @timed local_estimates(filtered, host_probes, data.bonds)
        density, bond_order = measurement.value.density, measurement.value.bond_order
        density_max, density_rms = metrics(density, source_density)
        bond_max, bond_rms = metrics(bond_order, source_bonds)
        result = (;
            density, bond_order, trace=sum(density),
            order=checkerboard_order(density, data, params.L),
            energy=hf_energy(data, density, bond_order, Float64(params.U)),
        )
        results[configuration.label] = result
        write_csv_row(io, (
            configuration.label, configuration.moments, configuration.probes,
            configuration.seed, result.trace, abs(result.trace - Ne), result.order,
            result.energy, kpm.time, measurement.time,
            density_max, density_rms, bond_max, bond_rms,
        ))
        filtered = nothing; probes_backend = nothing; host_probes = nothing
        GC.gc(); backend == :cuda && CUDA.reclaim()
    end
end

open(joinpath(output, "comparisons.csv"), "w") do io
    write_csv_row(io, (
        "left", "right", "interpretation", "density_max_abs_difference",
        "density_rms_difference", "bond_max_abs_difference", "bond_rms_difference",
        "checkerboard_order_difference", "energy_total_difference",
    ))
    comparisons = (
        ("M3200_R2048", "M3200_R4096", "nested_probe_count_at_fixed_M"),
        ("M3200_R2048", "M4000_R2048", "polynomial_order_at_fixed_nested_probes"),
        ("M4000_R2048", "M4000_R4096", "nested_probe_count_at_high_M"),
        ("M3200_R4096", "M4000_R4096", "polynomial_order_at_high_nested_probe_count"),
    )
    for (left_label, right_label, interpretation) in comparisons
        left, right = results[left_label], results[right_label]
        density_max, density_rms = metrics(left.density, right.density)
        bond_max, bond_rms = metrics(left.bond_order, right.bond_order)
        write_csv_row(io, (
            left_label, right_label, interpretation, density_max, density_rms,
            bond_max, bond_rms, left.order - right.order, left.energy - right.energy,
        ))
    end
end

println("Fixed-field KPM diagnostic complete: $output")
