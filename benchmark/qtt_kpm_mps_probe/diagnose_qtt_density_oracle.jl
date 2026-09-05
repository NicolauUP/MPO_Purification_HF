#!/usr/bin/env julia

"""Construct a QTT density directly from propagated MPS Hadamard probes.

The TCI oracle below never builds a density array. For an index `i` it only
contracts the already propagated MPS probes with `|i⟩`, then averages
`rα(i)⟨i|yα⟩`. A full array is materialised only *after* TCI for the small
fixed-H validation comparison.
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using TOML
using MPO_MeanField
using ITensors, ITensorMPS

length(ARGS) == 8 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS PROBES MAXDIM MPO_CUTOFF TCI_TOL",
)
campaign_file = abspath(ARGS[1]); task = parse(Int, ARGS[2]); output = abspath(ARGS[3])
moments = parse(Int, ARGS[4]); probes = parse(Int, ARGS[5]); maxdim = parse(Int, ARGS[6]); cutoff = parse(Float64, ARGS[7]); tci_tol = parse(Float64, ARGS[8])
isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments > 0 && probes > 0 && maxdim > 0 && cutoff > 0 && tci_tol > 0 || error("moments, probes, maxdim, cutoff, and tci_tol must be positive")

Base.include(Main, campaign_file)
@isdefined(campaign) || error("campaign file did not define `campaign`")
case = campaign_case(campaign, task); model = case.model
model isa SquareModel || error("diagnostic supports SquareModel only")
nx, ny = model.size; nx == ny || error("legacy QTT MPO construction supports equal squares only")
N = nx * ny; levels = qtt_levels(model); probes <= N || error("probes exceed N")

representation = QTTSettings(encoding=:interleaved, tci_tol=case.representation.tci_tol, cutoff=cutoff, maxdim=maxdim)
params = legacy_parameters(model, representation, case.solver)
system = System(params)
H_eff = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)

# Independent dense/vector reference, reordered to the QTT basis.
data = MPO_MeanField._kpm_data(model)
qtt_to_row = [begin x, y = square_lattice_decoder(index, levels); x + nx * y + 1 end for index in 0:(N - 1)]
H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
H = Matrix(H_row[qtt_to_row, qtt_to_row])
spectrum = eigen(Symmetric(H)); Ne = round(Int, model.filling * N)
center = (spectrum.values[1] + spectrum.values[end]) / 2
scale = (spectrum.values[end] - spectrum.values[1]) / 2 * 1.05
scaled_mu = ((spectrum.values[Ne] + spectrum.values[Ne + 1]) / 2 - center) / scale
coefficients = MPO_MeanField._kpm_coefficients(moments, scaled_mu)
H_scaled = +(
    H_eff / scale, (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=maxdim,
)

propagated = MPS[]
propagation_time = @elapsed for row in 0:(probes - 1)
    probe = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row)
    result = MPO_MeanField._qtt_mps_chebyshev_apply(H_scaled, probe, coefficients; cutoff=cutoff, maxdim=maxdim)
    push!(propagated, result.state)
end
probe_signs = [MPO_MeanField._qtt_hadamard_probe_vector(levels, row) for row in 0:(probes - 1)]
query_count = Ref(0)
function density_oracle(index)
    qtt_index = Int(index)
    query_count[] += 1
    return sum(probe_signs[column][qtt_index + 1] *
               MPO_MeanField._qtt_mps_amplitude(propagated[column], system.sites, qtt_index)
               for column in eachindex(propagated)) / probes
end

tci_time = @elapsed _, _, density_mps = MPO_MeanField.Quantics_TCI(density_oracle, Float64, system.sites, tci_tol)
tci_query_count = query_count[]

# Validation only: these arrays are deliberately materialised after fitting.
density_tci = MPO_MeanField._qtt_mps_amplitudes(density_mps, system.sites)
density_direct = [density_oracle(index) for index in 0:(N - 1)]
Z = hcat(probe_signs...)
PZ = MPO_MeanField._kpm_apply((H - center * Matrix{Float64}(I, N, N)) / scale, Z, coefficients)
density_vector = vec(sum(PZ .* Z; dims=2)) ./ probes

mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "matrix_dimension" => N, "moments" => moments, "probes" => probes,
        "maxdim" => maxdim, "mpo_cutoff" => cutoff, "tci_tolerance" => tci_tol, "propagation_time_s" => propagation_time,
        "tci_time_s" => tci_time, "tci_oracle_evaluations" => tci_query_count,
        "validation_oracle_evaluations" => query_count[] - tci_query_count,
        "density_tci_max_chi" => maxlinkdim(density_mps),
        "tci_vs_mps_oracle_max_abs_error" => maximum(abs.(density_tci - density_direct)),
        "tci_vs_mps_oracle_rms_error" => norm(density_tci - density_direct) / sqrt(N),
        "mps_oracle_vs_vector_max_abs_error" => maximum(abs.(density_direct - density_vector)),
        "mps_oracle_vs_vector_rms_error" => norm(density_direct - density_vector) / sqrt(N),
    ))
end
println("QTT density oracle diagnostic complete: $output")
