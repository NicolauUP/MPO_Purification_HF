#!/usr/bin/env julia

"""Benchmark a binary-QTT Hadamard probe register on a fixed 2D Hamiltonian.

For each `R = 2^q`, a q-bit probe register is interleaved with the q least
significant spatial QTT bits. This represents all Walsh columns
`z_a(i) = (-1)^popcount(i & a)` in one MPS, while the spatial MPO is extended
only by identity tensors on the register bits. The propagated register is
checked at sampled `(site, probe)` entries against the uncompressed sparse
vector KPM recurrence with exactly the same polynomial and spectral scaling.

This is a fixed-field representation benchmark, not a self-consistent HF run.
"""

using Dates
using LinearAlgebra
using Printf
using SparseArrays
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 8 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS R_VALUES MAXDIM CUTOFF BACKEND",
)
campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
probe_counts = parse.(Int, split(ARGS[5], ','))
maxdim = parse(Int, ARGS[6])
cutoff = parse(Float64, ARGS[7])
backend = Symbol(ARGS[8])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 4 || error("moments must be at least four")
!isempty(probe_counts) && issorted(probe_counts) && all(ispow2, probe_counts) ||
    error("R_VALUES must be a nonempty sorted list of powers of two")
maxdim > 0 || error("MAXDIM must be positive")
cutoff > 0 || error("CUTOFF must be positive")
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
    cutoff=cutoff, maxdim=maxdim,
)
parameters = legacy_parameters(model, representation, case.solver)
system = System(parameters)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)

# Use the same initial one-body Hamiltonian in the QTT MPO and the independent
# sparse-vector reference. Exact bounds are allowed for this N=4096 diagnostic.
data = MPO_MeanField._kpm_data(model)
H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
qtt_to_row = [begin
    x, y = square_lattice_decoder(index, levels)
    x + nx * y + 1
end for index in 0:(N - 1)]
H_qtt = sparse(H_row[qtt_to_row, qtt_to_row])
eigenvalues = eigvals(Symmetric(Matrix(H_qtt)))
Ne = round(Int, model.filling * N)
center = (first(eigenvalues) + last(eigenvalues)) / 2
scale = 1.05 * (last(eigenvalues) - first(eigenvalues)) / 2
mu = (eigenvalues[Ne] + eigenvalues[Ne + 1]) / 2
coefficients = MPO_MeanField._kpm_coefficients(moments, (mu - center) / scale)
H_scaled_vector = (H_qtt - center * sparse(I, N, N)) / scale
H_scaled_mpo = +(
    H_effective / scale,
    (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff, maxdim=maxdim,
)

mkpath(output)
progress_path = joinpath(output, "progress.log")
function progress(message)
    text = "[$(Dates.now())] $message"
    println(text); flush(stdout)
    open(progress_path, "a") do io; println(io, text); end
end
rms(values) = norm(values) / sqrt(length(values))
max_tensor_elements(psi::MPS) = maximum(tensor -> prod(dim(index) for index in inds(tensor)), psi)

# Deterministic spread across codes and spatial coordinates; this probes all
# register bits and both bulk and edge sites without an O(NR) host readback.
validation_indices = unique(round.(Int, range(0, N - 1; length=min(256, N))))
function validate_register(psi, register_sites, R, PZ)
    values = Matrix{Float64}(undef, length(validation_indices), R)
    for (row, index) in enumerate(validation_indices), code in 0:(R - 1)
        values[row, code + 1] = MPO_MeanField._qtt_probe_register_amplitude(
            psi, system.sites, register_sites, index, code,
        )
    end
    reference = PZ[validation_indices .+ 1, 1:R]
    signs = [isodd(count_ones(index & code)) ? -1.0 : 1.0
             for index in validation_indices, code in 0:(R - 1)]
    density_error = vec(sum((values - reference) .* signs; dims=2)) ./ R
    values - reference, density_error
end

# Compile the register-specific CUDA path outside the measured trajectories.
warm_register, warm_sites, _ = MPO_MeanField._qtt_hadamard_probe_register_mps(system.sites, 1)
warm_H = MPO_MeanField._qtt_extend_mpo_with_probe_register(H_scaled_mpo, system.sites, warm_sites)
warm_result = MPO_MeanField._qtt_mps_chebyshev_apply(
    ITensors.adapt(CUDA.CuArray, warm_H), ITensors.adapt(CUDA.CuArray, warm_register), coefficients[1:4];
    cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize(); warm_result = nothing; warm_H = nothing; warm_register = nothing
GC.gc(); CUDA.reclaim()

summary_rows = NamedTuple[]
for R in probe_counts
    probe_bits = trailing_zeros(R)
    progress("vector reference: R=$R, M=$moments")
    Z = [isodd(count_ones(index & code)) ? -1.0 : 1.0
         for index in 0:(N - 1), code in 0:(R - 1)]
    vector_time = @elapsed PZ = MPO_MeanField._kpm_apply(H_scaled_vector, Z, coefficients)

    register_host, register_sites, _ = MPO_MeanField._qtt_hadamard_probe_register_mps(
        system.sites, probe_bits,
    )
    extended_host = MPO_MeanField._qtt_extend_mpo_with_probe_register(
        H_scaled_mpo, system.sites, register_sites,
    )
    input_max_chi = maxlinkdim(register_host)
    input_max_elements = max_tensor_elements(register_host)
    progress("register CUDA propagation: R=$R (q=$probe_bits), χ=$maxdim")
    H_device = ITensors.adapt(CUDA.CuArray, extended_host)
    register_device = ITensors.adapt(CUDA.CuArray, register_host)
    CUDA.synchronize()
    gpu_time = @elapsed result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
        H_device, register_device, coefficients; cutoff=cutoff, maxdim=maxdim,
    )
    CUDA.synchronize()
    result_host = ITensors.cpu(result_device.state)
    difference, density_difference = validate_register(result_host, register_sites, R, PZ)
    push!(summary_rows, (
        probes=R, probe_bits=probe_bits, input_max_chi=input_max_chi,
        output_max_chi=maxlinkdim(result_host),
        output_mean_chi=MPO_MeanField._qtt_mps_mean_linkdim(result_host),
        input_max_tensor_elements=input_max_elements,
        output_max_tensor_elements=max_tensor_elements(result_host),
        vector_reference_time_s=vector_time, register_gpu_time_s=gpu_time,
        state_max_abs_difference=maximum(abs.(difference)),
        state_rms_difference=rms(difference),
        density_max_abs_difference=maximum(abs.(density_difference)),
        density_rms_difference=rms(density_difference),
    ))
    progress(@sprintf("register R=%d complete: %.3f s, output χ=%d, state RMS=%.3e, density RMS=%.3e",
        R, gpu_time, summary_rows[end].output_max_chi,
        summary_rows[end].state_rms_difference, summary_rows[end].density_rms_difference))
    H_device = nothing; register_device = nothing; result_device = nothing
    extended_host = nothing; register_host = nothing; result_host = nothing; PZ = nothing; Z = nothing
    GC.gc(); CUDA.reclaim()
end

open(joinpath(output, "summary.csv"), "w") do io
    println(io, "probes,probe_bits,input_max_chi,output_max_chi,output_mean_chi,input_max_tensor_elements,output_max_tensor_elements,vector_reference_time_s,register_gpu_time_s,state_max_abs_difference,state_rms_difference,density_max_abs_difference,density_rms_difference")
    for row in summary_rows
        println(io, join(values(row), ','))
    end
end
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file, "label" => case.label,
        "diagnostic" => "qtt_binary_probe_register_vs_sparse_vector_kpm",
        "matrix_dimension" => N, "moments" => moments, "probe_counts" => probe_counts,
        "maxdim" => maxdim, "cutoff" => cutoff, "backend" => "cuda",
        "device_name" => CUDA.name(CUDA.device()), "validation_sites" => length(validation_indices),
        "probe_family" => "standard Walsh codes (-1)^popcount(index & code)",
        "comparison" => "register amplitudes versus uncompressed sparse-vector KPM",
    ))
end
progress("binary probe-register benchmark complete: $output")
