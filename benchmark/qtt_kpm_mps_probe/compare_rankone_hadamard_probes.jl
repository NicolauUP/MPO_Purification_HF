#!/usr/bin/env julia

"""Compare rank-one QTT Hadamard probes with production randomized probes.

This fixed-H diagnostic intentionally uses the same Jackson polynomial and
exact dense reference for both probe families. It answers one question only:
does removing the random diagonal sign field (so every probe is rank-one in
QTT) materially worsen finite-R local density/bond estimates?
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using TOML
using MPO_MeanField

length(ARGS) == 7 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS PROBE_LIST SEED BACKEND",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
probe_counts = parse.(Int, split(ARGS[5], ','))
seed = parse(Int, ARGS[6])
backend = Symbol(lowercase(ARGS[7]))
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
backend in (:cpu, :cuda) || error("backend must be cpu or cuda")
all(>(0), probe_counts) || error("every probe count must be positive")

# Include before any model-dependent work so campaign closures are valid in
# TCI/MPO construction and in the KPM data constructor.
Base.include(Main, campaign_file)
@isdefined(campaign) || error("campaign file did not define `campaign`")
case = campaign_case(campaign, task)
model = case.model
model isa SquareModel || error("diagnostic supports SquareModel only")
nx, ny = model.size
nx == ny || error("rank-one QTT comparison currently supports equal squares only")
N = nx * ny
maximum(probe_counts) <= N || error("probe count exceeds N=$N")
levels = qtt_levels(model)

qtt_to_row = [begin
    x, y = square_lattice_decoder(index, levels)
    x + nx * y + 1
end for index in 0:(N - 1)]
row_to_qtt = invperm(qtt_to_row)
data = MPO_MeanField._kpm_data(model)
H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
H = Matrix(H_row[qtt_to_row, qtt_to_row])
bonds = [(row_to_qtt[i], row_to_qtt[j]) for (i, j, _) in data.bonds]
Ne = round(Int, model.filling * N)
spectrum = eigen(Symmetric(H))
μ = (spectrum.values[Ne] + spectrum.values[Ne + 1]) / 2
center = (spectrum.values[1] + spectrum.values[end]) / 2
scale = (spectrum.values[end] - spectrum.values[1]) / 2 * 1.05
coefficients = MPO_MeanField._kpm_coefficients(moments, (μ - center) / scale)
scaled_H = (H - center * Matrix{Float64}(I, N, N)) / scale
occupied = spectrum.vectors[:, 1:Ne]
exact_density = vec(sum(abs2, occupied; dims=2))
exact_bonds = [dot(@view(occupied[i, :]), @view(occupied[j, :])) for (i, j) in bonds]

cuda = nothing
if backend == :cuda
    @eval using CUDA
    CUDA.functional() || error("CUDA is not functional on this node")
    cuda = CUDA
    scaled_H_backend = CUDA.CuArray(scaled_H)
else
    scaled_H_backend = scaled_H
end

function local_observables(PZ, Z, bonds)
    R = size(Z, 2)
    density = vec(sum(PZ .* Z; dims=2)) ./ R
    order = [
        (dot(@view(PZ[i, :]), @view(Z[j, :])) + dot(@view(PZ[j, :]), @view(Z[i, :]))) / (2R)
        for (i, j) in bonds
    ]
    trace = sum(PZ .* Z) / R
    return (; density, order, trace)
end

function relative_rms(value, reference)
    return norm(value - reference) / sqrt(length(value))
end

function run_family(name, probes)
    device_probes = isnothing(cuda) ? probes : cuda.CuArray(probes)
    elapsed = @elapsed filtered_device = MPO_MeanField._kpm_apply(
        scaled_H_backend, device_probes, coefficients;
        synchronize=isnothing(cuda) ? (() -> nothing) : cuda.synchronize,
    )
    filtered = isnothing(cuda) ? filtered_device : Array(filtered_device)
    estimates = local_observables(filtered, probes, bonds)
    return (; name, elapsed, estimates)
end

mkpath(output)
rows = Vector{Vector{String}}()
for R in probe_counts
    rankone = hcat([
        MPO_MeanField._qtt_hadamard_probe_vector(levels, row)
        for row in 0:(R - 1)
    ]...)
    # Production probes are generated in row-major physical storage and then
    # permuted into the QTT basis used by the rank-one rows above.
    randomized = MPO_MeanField._kpm_hadamard(data.codes, R, seed)[qtt_to_row, :]
    deterministic_result = run_family("rankone", rankone)
    randomized_result = run_family("randomized", randomized)
    for result in (deterministic_result, randomized_result)
        push!(rows, [
            string(R), result.name, @sprintf("%.16g", result.elapsed),
            @sprintf("%.16g", result.estimates.trace),
            @sprintf("%.16g", maximum(abs.(result.estimates.density - exact_density))),
            @sprintf("%.16g", relative_rms(result.estimates.density, exact_density)),
            @sprintf("%.16g", maximum(abs.(result.estimates.order - exact_bonds))),
            @sprintf("%.16g", relative_rms(result.estimates.order, exact_bonds)),
        ])
    end
end
open(joinpath(output, "probe_comparison.csv"), "w") do io
    println(io, "probes,family,kpm_time_s,trace,density_max_abs_error,density_rms_error,bond_max_abs_error,bond_rms_error")
    for row in rows
        println(io, join(row, ','))
    end
end
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "matrix_dimension" => N, "moments" => moments, "probe_counts" => probe_counts,
        "seed" => seed, "backend" => string(backend), "spectral_lower_exact" => spectrum.values[1],
        "spectral_upper_exact" => spectrum.values[end], "exact_fermi_gap" => spectrum.values[Ne + 1] - spectrum.values[Ne],
        "probe_families" => ["rankone_qtt_hadamard", "randomized_production_hadamard"],
    ))
end
println("rank-one versus randomized Hadamard comparison complete: $output")
