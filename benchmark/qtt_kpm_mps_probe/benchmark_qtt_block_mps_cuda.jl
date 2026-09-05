#!/usr/bin/env julia

"""Benchmark batched (block-MPS) versus independent QTT KPM probe propagation.

The block state has one spectator index labelling Hadamard probes. The
Hamiltonian MPO acts only on quantics sites, so each block column represents
the same Chebyshev recurrence as a separately propagated probe. We compare
the resulting amplitudes to the independent GPU trajectories, while reporting
the extra rank/cap required by the shared block representation.
"""

using Dates
using LinearAlgebra
using Printf
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 9 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS R_VALUES SEQ_MAXDIM BLOCK_MAXDIMS CUTOFF BACKEND",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
probe_counts = parse.(Int, split(ARGS[5], ','))
sequential_maxdim = parse(Int, ARGS[6])
block_maxdims = parse.(Int, split(ARGS[7], ','))
cutoff = parse(Float64, ARGS[8])
backend = Symbol(ARGS[9])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 2 || error("moments must be at least two")
!isempty(probe_counts) && issorted(probe_counts) && all(>(0), probe_counts) ||
    error("R_VALUES must be a nonempty sorted list of positive integers")
sequential_maxdim > 0 || error("SEQ_MAXDIM must be positive")
!isempty(block_maxdims) && all(>(0), block_maxdims) || error("BLOCK_MAXDIMS must be positive")
cutoff > 0 || error("cutoff must be positive")
backend == :cuda || error("this benchmark currently requires backend=cuda")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
case = campaign_case(Main.campaign, task)
model = case.model
model isa SquareModel || error("diagnostic requires SquareModel")
nx, ny = model.size
nx == ny || error("diagnostic currently requires a square model")
N = nx * ny
levels = qtt_levels(model)
Rmax = maximum(probe_counts)
Rmax <= N || error("maximum probe count $Rmax exceeds N=$N")

representation = QTTSettings(
    encoding=:interleaved, tci_tol=case.representation.tci_tol,
    cutoff=cutoff, maxdim=maximum((sequential_maxdim, block_maxdims...)),
)
parameters = legacy_parameters(model, representation, case.solver)
system = System(parameters)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=representation.maxdim)

# Exact bounds are permitted here because this is a controlled N=1024
# diagnostic; no spectrum-estimation issue is being benchmarked.
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
scale = 1.05 * (spectrum.values[end] - spectrum.values[1]) / 2
mu = (spectrum.values[Ne] + spectrum.values[Ne + 1]) / 2
coefficients = MPO_MeanField._kpm_coefficients(moments, (mu - center) / scale)
H_scaled_host = +(
    H_effective / scale,
    (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=representation.maxdim,
)
H_scaled_device = ITensors.adapt(CUDA.CuArray, H_scaled_host)
rows = collect(0:(Rmax - 1))

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    text = "[$(Dates.now())] $message"
    println(text); flush(stdout)
    open(progress_path, "a") do io; println(io, text); end
end
rms(values) = norm(values) / sqrt(length(values))
function max_tensor_elements(psi::MPS)
    maximum(tensor -> prod(dim(index) for index in inds(tensor)), psi)
end

# Warm compilation off the timed path. This uses a degree-two polynomial and
# a single rank-one column, not one of the measured trajectories.
warm_probe = ITensors.adapt(CUDA.CuArray, MPO_MeanField._qtt_hadamard_probe_mps(system.sites, 0))
MPO_MeanField._qtt_mps_chebyshev_apply(H_scaled_device, warm_probe, coefficients[1:2];
    cutoff=cutoff, maxdim=sequential_maxdim)
CUDA.synchronize(); warm_probe = nothing; CUDA.reclaim()

progress("sequential GPU baseline: N=$N M=$moments Rmax=$Rmax χ=$sequential_maxdim")
sequential_states = MPS[]
sequential_time = 0.0
sequential_rows = NamedTuple[]
for (ordinal, row) in enumerate(rows)
    probe_device = ITensors.adapt(CUDA.CuArray, MPO_MeanField._qtt_hadamard_probe_mps(system.sites, row))
    CUDA.synchronize()
    elapsed = @elapsed result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_scaled_device, probe_device, coefficients; cutoff=cutoff, maxdim=sequential_maxdim,
    )
    CUDA.synchronize()
    state_host = ITensors.cpu(result_device.state)
    push!(sequential_states, state_host)
    push!(sequential_rows, (ordinal=ordinal, row=row, time_s=elapsed,
        max_chi=maxlinkdim(state_host), mean_chi=MPO_MeanField._qtt_mps_mean_linkdim(state_host),
        max_tensor_elements=max_tensor_elements(state_host)))
    # This diagnostic is intentionally a top-level script. Julia's global
    # loop scope therefore requires an explicit global update for the
    # accumulated baseline time.
    global sequential_time += elapsed
    probe_device = nothing; result_device = nothing
    GC.gc(); CUDA.reclaim()
    progress(@sprintf("sequential %d/%d: %.3f s, χ=%d", ordinal, Rmax, elapsed, maxlinkdim(state_host)))
end

# Full amplitude extraction is cheap at N=1024. For the 64x64 scaling test,
# compare a fixed spread of sites instead: its purpose is to diagnose the
# block representation, not to pay an O(NR) host contraction cost.
validation_indices = N <= 1024 ? collect(0:(N - 1)) : unique(round.(Int,
    range(0, N - 1; length=min(512, N))))
sequential_values = hcat([[
    MPO_MeanField._qtt_mps_amplitude(state, system.sites, index)
    for index in validation_indices
] for state in sequential_states]...)
Z = hcat([MPO_MeanField._qtt_hadamard_probe_vector(levels, row)[validation_indices .+ 1]
    for row in rows]...)

summary_rows = NamedTuple[]
for R in probe_counts, block_maxdim in block_maxdims
    active_rows = rows[1:R]
    progress("block GPU propagation: R=$R χ=$block_maxdim")
    block_host, probe_index = MPO_MeanField._qtt_hadamard_probe_block_mps(
        system.sites, active_rows; cutoff=cutoff, maxdim=block_maxdim,
    )
    input_max_chi = maxlinkdim(block_host)
    input_max_elements = max_tensor_elements(block_host)
    block_device = ITensors.adapt(CUDA.CuArray, block_host)
    CUDA.synchronize()
    elapsed = @elapsed result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_scaled_device, block_device, coefficients; cutoff=cutoff, maxdim=block_maxdim,
    )
    CUDA.synchronize()
    result_host = ITensors.cpu(result_device.state)
    block_values = hcat([[
        MPO_MeanField._qtt_block_mps_amplitude(result_host, system.sites, probe_index, slot, index)
        for index in validation_indices
    ] for slot in 1:R]...)
    block_density = vec(sum(block_values .* Z[:, 1:R]; dims=2)) ./ R
    sequential_density = vec(sum(sequential_values[:, 1:R] .* Z[:, 1:R]; dims=2)) ./ R
    difference = block_values - sequential_values[:, 1:R]
    density_difference = block_density - sequential_density
    push!(summary_rows, (
        probes=R, block_maxdim=block_maxdim, input_max_chi=input_max_chi,
        output_max_chi=maxlinkdim(result_host),
        output_mean_chi=MPO_MeanField._qtt_mps_mean_linkdim(result_host),
        input_max_tensor_elements=input_max_elements,
        output_max_tensor_elements=max_tensor_elements(result_host),
        block_time_s=elapsed, sequential_total_time_s=sequential_time,
        speedup=sequential_time / elapsed,
        state_max_abs_difference=maximum(abs.(difference)),
        state_rms_difference=rms(difference),
        density_max_abs_difference=maximum(abs.(density_difference)),
        density_rms_difference=rms(density_difference),
    ))
    block_device = nothing; result_device = nothing; result_host = nothing
    GC.gc(); CUDA.reclaim()
    progress(@sprintf("block R=%d χ=%d complete: %.3f s, output χ=%d, state RMS=%.3e, density RMS=%.3e",
        R, block_maxdim, elapsed, summary_rows[end].output_max_chi,
        summary_rows[end].state_rms_difference, summary_rows[end].density_rms_difference))
end

open(joinpath(output, "sequential.csv"), "w") do io
    println(io, "ordinal,probe_row,time_s,max_chi,mean_chi,max_tensor_elements")
    for row in sequential_rows
        println(io, "$(row.ordinal),$(row.row),$(row.time_s),$(row.max_chi),$(row.mean_chi),$(row.max_tensor_elements)")
    end
end
open(joinpath(output, "summary.csv"), "w") do io
    println(io, "probes,block_maxdim,input_max_chi,output_max_chi,output_mean_chi,input_max_tensor_elements,output_max_tensor_elements,block_time_s,sequential_prefix_time_s,sequential_over_block_speedup,state_max_abs_difference,state_rms_difference,density_max_abs_difference,density_rms_difference")
    for row in summary_rows
        prefix_time = sum(item.time_s for item in sequential_rows[1:row.probes])
        println(io, "$(row.probes),$(row.block_maxdim),$(row.input_max_chi),$(row.output_max_chi),$(row.output_mean_chi),$(row.input_max_tensor_elements),$(row.output_max_tensor_elements),$(row.block_time_s),$(prefix_time),$(prefix_time / row.block_time_s),$(row.state_max_abs_difference),$(row.state_rms_difference),$(row.density_max_abs_difference),$(row.density_rms_difference)")
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "diagnostic" => "qtt_block_mps_vs_independent_gpu_kpm",
        "matrix_dimension" => N, "moments" => moments, "probe_counts" => probe_counts,
        "sequential_maxdim" => sequential_maxdim, "block_maxdims" => block_maxdims,
        "cutoff" => cutoff, "backend" => "cuda", "device_name" => CUDA.name(CUDA.device()),
        "comparison" => "block columns versus independently propagated columns; matching cutoff",
        "validation_sites" => length(validation_indices),
    ))
end
progress("block-MPS benchmark complete: $output")
