#!/usr/bin/env julia

using LinearAlgebra
using TOML

length(ARGS) == 3 || error(
    "usage: compare_kpm_dense_outputs.jl KPM_DIRECTORY DENSE_DIRECTORY OUTPUT",
)
kpm_directory = abspath(ARGS[1])
dense_directory = abspath(ARGS[2])
output = abspath(ARGS[3])

parse_fields(line) = split(replace(chomp(line), "\"" => ""), ',')

function load_keyed_column(path, key_columns, value_column)
    values = Dict{Tuple,Float64}()
    open(path) do io
        header = parse_fields(readline(io))
        key_positions = map(key_columns) do column
            position = findfirst(==(column), header)
            isnothing(position) &&
                error("missing key column $column in $path")
            position
        end
        value_position = findfirst(==(value_column), header)
        isnothing(value_position) &&
            error("missing value column $value_column in $path")
        for line in eachline(io)
            fields = parse_fields(line)
            key = Tuple(fields[position] for position in key_positions)
            haskey(values, key) && error("duplicate key $key in $path")
            values[key] = parse(Float64, fields[value_position])
        end
    end
    return values
end

kpm_density = load_keyed_column(
    joinpath(kpm_directory, "density.csv"), ("site",), "audit",
)
dense_density = load_keyed_column(
    joinpath(dense_directory, "site_density.csv"), ("site",), "density",
)
kpm_bonds = load_keyed_column(
    joinpath(kpm_directory, "bond_order.csv"),
    ("site", "neighbour", "orientation"),
    "audit",
)
dense_bonds = load_keyed_column(
    joinpath(dense_directory, "bond_order.csv"),
    ("site_left", "site_right", "orientation"),
    "real",
)
Set(keys(kpm_density)) == Set(keys(dense_density)) ||
    error("density keys differ")
Set(keys(kpm_bonds)) == Set(keys(dense_bonds)) ||
    error("bond keys differ")

kpm_observables = TOML.parsefile(joinpath(kpm_directory, "observables.toml"))
dense_observables =
    TOML.parsefile(joinpath(dense_directory, "observables.toml"))
dense_metadata =
    TOML.parsefile(joinpath(dense_directory, "metadata.toml"))
density_difference = [
    kpm_density[key] - dense_density[key] for key in keys(kpm_density)
]
bond_difference = [
    kpm_bonds[key] - dense_bonds[key] for key in keys(kpm_bonds)
]

summary = Dict(
    "matrix_dimension" => length(kpm_density),
    "density_max_abs_error" => maximum(abs, density_difference),
    "density_rms_error" =>
        norm(density_difference) / sqrt(length(density_difference)),
    "bond_max_abs_error" => maximum(abs, bond_difference),
    "bond_rms_error" =>
        norm(bond_difference) / sqrt(length(bond_difference)),
    "particle_number_difference" =>
        kpm_observables["particle_number"] -
        dense_observables["particle_number"],
    "energy_total_difference" =>
        kpm_observables["audited_energy_total"] -
        dense_observables["energy_total"],
    "kpm_scf_converged" => kpm_observables["scf_converged"],
    "dense_scf_converged" => dense_metadata["scf_converged"],
    "kpm_directory" => kpm_directory,
    "dense_directory" => dense_directory,
)
open(output, "w") do io
    TOML.print(io, summary)
end
println("KPM/dense comparison written to $output")
