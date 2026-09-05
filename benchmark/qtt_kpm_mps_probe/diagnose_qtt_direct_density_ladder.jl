#!/usr/bin/env julia

"""Validate direct QTT density-MPS assembly for a nested Hadamard probe ladder.

The MPS density is assembled as mean(zα .* Pzα), not fitted through TCI.
An independent vector KPM recurrence and dense fixed-H projector at 32x32
separate assembly error, MPO--MPS propagation error, and finite-R/KPM error.
"""

using Dates
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 7 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS MAXDIM CUTOFF R_VALUES",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
maxdim = parse(Int, ARGS[5])
cutoff = parse(Float64, ARGS[6])
probe_counts = parse.(Int, split(ARGS[7], ','))

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 2 || error("moments must be at least two")
maxdim > 0 || error("maxdim must be positive")
cutoff > 0 || error("cutoff must be positive")
!isempty(probe_counts) && issorted(probe_counts) && all(>(0), probe_counts) ||
    error("R_VALUES must be positive and sorted")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
case = campaign_case(Main.campaign, task)
model = case.model
model isa SquareModel || error("diagnostic requires SquareModel")
nx, ny = model.size
nx == ny || error("diagnostic currently requires an equal square")
N = nx * ny
levels = qtt_levels(model)
maximum(probe_counts) <= N || error("probe count exceeds N=$N")

representation = QTTSettings(
    encoding=:interleaved, tci_tol=case.representation.tci_tol,
    cutoff=cutoff, maxdim=maxdim,
)
parameters = legacy_parameters(model, representation, case.solver)
system = System(parameters)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)

data = MPO_MeanField._kpm_data(model)
H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
qtt_to_row = [begin
    x, y = square_lattice_decoder(index, levels)
    x + nx * y + 1
end for index in 0:(N - 1)]
H_qtt = Matrix(H_row[qtt_to_row, qtt_to_row])
spectrum = eigen(Symmetric(H_qtt))
Ne = round(Int, model.filling * N)
center = (spectrum.values[1] + spectrum.values[end]) / 2
scale = (spectrum.values[end] - spectrum.values[1]) / 2 * 1.05
coefficients = MPO_MeanField._kpm_coefficients(moments,
    ((spectrum.values[Ne] + spectrum.values[Ne + 1]) / 2 - center) / scale)
H_scaled_host = +(
    H_effective / scale,
    (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=maxdim,
)
H_scaled_device = ITensors.adapt(CUDA.CuArray, H_scaled_host)
H_scaled_dense = (H_qtt - center * Matrix{Float64}(I, N, N)) / scale

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    text = "[$(Dates.now())] $message"
    println(text); flush(stdout)
    open(progress_path, "a") do io; println(io, text); end
end

Rmax = maximum(probe_counts)
rows = collect(0:(Rmax - 1))
states = MPS[]
propagation_rows = NamedTuple[]
progress("starting CUDA MPS propagation: N=$N M=$moments Rmax=$Rmax")
for (ordinal, row) in enumerate(rows)
    probe_host = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row)
    probe_device = ITensors.adapt(CUDA.CuArray, probe_host)
    CUDA.synchronize()
    elapsed = @elapsed state_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_scaled_device, probe_device, coefficients; cutoff=cutoff, maxdim=maxdim,
    )
    CUDA.synchronize()
    state_host = ITensors.cpu(state_device.state)
    push!(states, state_host)
    final_chi = maxlinkdim(state_host)
    probe_device = nothing; state_device = nothing
    GC.gc(); CUDA.reclaim()
    push!(propagation_rows, (ordinal=ordinal, row=row, time_s=elapsed, max_chi=final_chi))
    progress(@sprintf("probe %d/%d complete: %.3f s, χ=%d", ordinal, Rmax, elapsed, final_chi))
end

# Independent vector KPM and exact dense projector on the same H and QTT order.
Z = hcat([MPO_MeanField._qtt_hadamard_probe_vector(levels, row) for row in rows]...)
PZ = MPO_MeanField._kpm_apply(H_scaled_dense, Z, coefficients)
exact_projector = spectrum.vectors[:, 1:Ne] * spectrum.vectors[:, 1:Ne]'
exact_density = real.(diag(exact_projector))
mps_outputs = hcat([MPO_MeanField._qtt_mps_amplitudes(state, system.sites) for state in states]...)

open(joinpath(output, "propagation.csv"), "w") do io
    println(io, "ordinal,probe_row,time_s,final_max_chi")
    for row in propagation_rows
        println(io, "$(row.ordinal),$(row.row),$(row.time_s),$(row.max_chi)")
    end
end

rms(value) = norm(value) / sqrt(length(value))
summary_rows = NamedTuple[]
for R in probe_counts
    density_time = @elapsed direct_density_mps = MPO_MeanField._qtt_density_mps_from_hadamard_probes(
        states[1:R], system.sites, rows[1:R]; cutoff=cutoff, maxdim=maxdim,
    )
    direct_density = MPO_MeanField._qtt_mps_amplitudes(direct_density_mps, system.sites)
    raw_mps_density = vec(sum(mps_outputs[:, 1:R] .* Z[:, 1:R]; dims=2)) ./ R
    vector_density = vec(sum(PZ[:, 1:R] .* Z[:, 1:R]; dims=2)) ./ R
    push!(summary_rows, (
        probes=R, density_max_chi=maxlinkdim(direct_density_mps), assembly_time_s=density_time,
        direct_trace=sum(direct_density),
        direct_vs_raw_mps_max=maximum(abs.(direct_density - raw_mps_density)),
        direct_vs_raw_mps_rms=rms(direct_density - raw_mps_density),
        direct_vs_vector_max=maximum(abs.(direct_density - vector_density)),
        direct_vs_vector_rms=rms(direct_density - vector_density),
        vector_vs_exact_max=maximum(abs.(vector_density - exact_density)),
        vector_vs_exact_rms=rms(vector_density - exact_density),
    ))
    progress(@sprintf("R=%d direct density complete: χ=%d, assembly RMS=%.3e, vector/ED RMS=%.3e",
        R, summary_rows[end].density_max_chi, summary_rows[end].direct_vs_raw_mps_rms,
        summary_rows[end].vector_vs_exact_rms))
end

open(joinpath(output, "summary.csv"), "w") do io
    println(io, "probes,density_mps_max_chi,assembly_time_s,direct_trace,direct_vs_raw_mps_max_abs_error,direct_vs_raw_mps_rms_error,direct_vs_vector_max_abs_error,direct_vs_vector_rms_error,vector_vs_exact_max_abs_error,vector_vs_exact_rms_error")
    for row in summary_rows
        println(io, "$(row.probes),$(row.density_max_chi),$(row.assembly_time_s),$(row.direct_trace),$(row.direct_vs_raw_mps_max),$(row.direct_vs_raw_mps_rms),$(row.direct_vs_vector_max),$(row.direct_vs_vector_rms),$(row.vector_vs_exact_max),$(row.vector_vs_exact_rms)")
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "matrix_dimension" => N, "moments" => moments, "probe_counts" => probe_counts,
        "maxdim" => maxdim, "cutoff" => cutoff, "device_name" => CUDA.name(CUDA.device()),
        "construction" => "direct finite-probe density MPS sum; no TCI",
    ))
end
progress("direct density MPS ladder complete: $output")
