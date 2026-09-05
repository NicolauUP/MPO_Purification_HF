#!/usr/bin/env julia

"""Large fixed-H direct QTT density/Hartree ladder with GPU MPS propagation.

For nested rank-one Hadamard probe counts, construct the density directly as
an MPS and obtain the Hartree field by applying the exact QTT binary-carry
nearest-neighbour adjacency MPO. No TCI fit and no full N-site density array
are constructed. Deterministic sample points validate both QTT fields against
the same uncompressed finite-probe MPS oracle.
"""

using Dates
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 8 || error(
    "usage: $(PROGRAM_FILE) LEGACY_CAMPAIGN.jl TASK OUTPUT MOMENTS MAXDIM CUTOFF R_VALUES SAMPLE_COUNT",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
maxdim = parse(Int, ARGS[5])
cutoff = parse(Float64, ARGS[6])
probe_counts = parse.(Int, split(ARGS[7], ','))
sample_count = parse(Int, ARGS[8])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 2 || error("moments must be at least two")
maxdim > 0 || error("maxdim must be positive")
cutoff > 0 || error("cutoff must be positive")
!isempty(probe_counts) && issorted(probe_counts) && all(>(0), probe_counts) ||
    error("R_VALUES must be positive and sorted")
sample_count > 0 || error("SAMPLE_COUNT must be positive")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
spec = Main.campaign.runs[task]
params = spec.params
params isa ParametersSquare || error("diagnostic requires ParametersSquare")
N = 1 << params.L
Rmax = maximum(probe_counts)
Rmax <= N || error("probe count exceeds N=$N")
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
abs(lower + upper) <= 1e-12 || error("large direct QTT diagnostic requires symmetric spectral bounds")

effective_params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=cutoff, itensors_maxdim=maxdim,
    density=params.density, purification_steps=params.purification_steps,
    scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
    scf_max_iterations=params.scf_max_iterations,
)
system = System(effective_params)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)
H_device = ITensors.adapt(CUDA.CuArray, H_effective / ((upper - lower) / 2))
coefficients = MPO_MeanField._kpm_coefficients(moments, 0.0)

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    text = "[$(Dates.now())] $message"
    println(text); flush(stdout)
    open(progress_path, "a") do io; println(io, text); end
end

rows = collect(0:(Rmax - 1))
states = MPS[]
propagation_rows = NamedTuple[]
progress("starting GPU propagation: N=$N M=$moments Rmax=$Rmax χmax=$maxdim")
for (ordinal, row) in enumerate(rows)
    probe_device = ITensors.adapt(CUDA.CuArray,
        MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row))
    CUDA.synchronize()
    elapsed = @elapsed state_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_device, probe_device, coefficients; cutoff=cutoff, maxdim=maxdim,
    )
    CUDA.synchronize()
    state_host = ITensors.cpu(state_device.state)
    push!(states, state_host)
    final_chi = maxlinkdim(state_host)
    probe_device = nothing; state_device = nothing
    GC.gc(); CUDA.reclaim()
    free_gpu = CUDA.free_memory()
    push!(propagation_rows, (ordinal=ordinal, row=row, time_s=elapsed,
        final_max_chi=final_chi, free_gpu=free_gpu))
    progress(@sprintf("probe %d/%d complete: %.3f s, χ=%d, free GPU %.2f GiB",
        ordinal, Rmax, elapsed, final_chi, free_gpu / 2.0^30))
end

signs = [MPO_MeanField._qtt_hadamard_probe_vector(params.L, row) for row in rows]
validation_indices = unique(round.(Int, range(0, N - 1; length=min(sample_count, N))))
ones_mps = MPS([ITensor([1.0, 1.0], site) for site in system.sites])
adjacency = square_neighbour_adjacency_mpo(system.sites)

function raw_density(index::Int, R::Int)
    return sum(signs[column][index + 1] *
               MPO_MeanField._qtt_mps_amplitude(states[column], system.sites, index)
               for column in 1:R) / R
end
function raw_hartree(index::Int, R::Int)
    site = index + 1
    return params.U * sum(raw_density(neighbour - 1, R)
        for neighbour in values(square_neighbours(site, params.L)) if !isnothing(neighbour))
end

open(joinpath(output, "propagation.csv"), "w") do io
    println(io, "ordinal,probe_row,time_s,final_max_chi,free_gpu_bytes")
    for row in propagation_rows
        println(io, "$(row.ordinal),$(row.row),$(row.time_s),$(row.final_max_chi),$(row.free_gpu)")
    end
end

summary_rows = NamedTuple[]
validation_data = Dict{Int,NamedTuple}()
for R in probe_counts
    density_time = @elapsed density = MPO_MeanField._qtt_density_mps_from_hadamard_probes(
        states[1:R], system.sites, rows[1:R]; cutoff=cutoff, maxdim=maxdim,
    )
    hartree_time = @elapsed hartree = params.U * apply(adjacency, density;
        cutoff=cutoff, maxdim=maxdim)
    raw_density_values = [raw_density(index, R) for index in validation_indices]
    raw_hartree_values = [raw_hartree(index, R) for index in validation_indices]
    density_values = [MPO_MeanField._qtt_mps_amplitude(density, system.sites, index)
                      for index in validation_indices]
    hartree_values = [MPO_MeanField._qtt_mps_amplitude(hartree, system.sites, index)
                      for index in validation_indices]
    density_error = density_values .- raw_density_values
    hartree_error = hartree_values .- raw_hartree_values
    row = (
        probes=R,
        density_max_chi=maxlinkdim(density),
        hartree_max_chi=maxlinkdim(hartree),
        density_time_s=density_time,
        hartree_time_s=hartree_time,
        density_trace=real(inner(ones_mps, density)),
        density_trace_error=abs(real(inner(ones_mps, density)) - params.density * N),
        density_sample_maximum=maximum(abs, density_error),
        density_sample_rms=norm(density_error) / sqrt(length(density_error)),
        hartree_sample_maximum=maximum(abs, hartree_error),
        hartree_sample_rms=norm(hartree_error) / sqrt(length(hartree_error)),
    )
    push!(summary_rows, row)
    validation_data[R] = (raw_density=raw_density_values, density=density_values,
        raw_hartree=raw_hartree_values, hartree=hartree_values)
    progress(@sprintf("R=%d complete: χ(n)=%d χ(VH)=%d, density sample RMS=%.3e, Hartree sample RMS=%.3e",
        R, row.density_max_chi, row.hartree_max_chi,
        row.density_sample_rms, row.hartree_sample_rms))
end

open(joinpath(output, "summary.csv"), "w") do io
    println(io, "probes,density_max_chi,hartree_max_chi,density_time_s,hartree_time_s,density_trace,density_trace_error,density_sample_max_abs_error,density_sample_rms_error,hartree_sample_max_abs_error,hartree_sample_rms_error")
    for row in summary_rows
        println(io, "$(row.probes),$(row.density_max_chi),$(row.hartree_max_chi),$(row.density_time_s),$(row.hartree_time_s),$(row.density_trace),$(row.density_trace_error),$(row.density_sample_maximum),$(row.density_sample_rms),$(row.hartree_sample_maximum),$(row.hartree_sample_rms)")
    end
end
open(joinpath(output, "validation.csv"), "w") do io
    header = ["qtt_index"]
    for R in probe_counts
        append!(header, ["raw_density_R$R", "density_R$R", "raw_hartree_R$R", "hartree_R$R"])
    end
    println(io, join(header, ','))
    for position in eachindex(validation_indices)
        row = [string(validation_indices[position])]
        for R in probe_counts
            values = validation_data[R]
            append!(row, [@sprintf("%.16g", values.raw_density[position]),
                @sprintf("%.16g", values.density[position]),
                @sprintf("%.16g", values.raw_hartree[position]),
                @sprintf("%.16g", values.hartree[position])])
        end
        println(io, join(row, ','))
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => spec.label,
        "matrix_dimension" => N, "moments" => moments, "probe_counts" => probe_counts,
        "maxdim" => maxdim, "cutoff" => cutoff, "sample_count" => length(validation_indices),
        "device_name" => CUDA.name(CUDA.device()), "spectral_lower" => lower,
        "spectral_upper" => upper,
        "construction" => "direct finite-probe density MPS followed by binary-carry adjacency MPO Hartree application; no TCI",
        "validation_note" => "sampled against the same uncompressed finite-probe MPS oracle; no ED/vector reference at 256x256",
    ))
end
progress("direct density/Hartree ladder complete: $output")
