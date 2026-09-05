#!/usr/bin/env julia

"""Compare one fixed-H QTT MPO--MPS KPM trajectory on CPU and CUDA.

This is a backend-validation benchmark, not an SCF calculation. It keeps the
model, spectral scaling, rank-one Hadamard probe, Chebyshev coefficients,
cutoff, and maximum bond dimension identical between the two paths. The CUDA
timing excludes one short warm-up trajectory used solely to compile kernels.
"""

using Dates
using LinearAlgebra
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 7 || error(
    "usage: $(PROGRAM_FILE) CAMPAIGN.jl TASK OUTPUT MOMENTS MAXDIM CUTOFF PROBE_ROW",
)

campaign_file = abspath(ARGS[1])
task = parse(Int, ARGS[2])
output = abspath(ARGS[3])
moments = parse(Int, ARGS[4])
maxdim = parse(Int, ARGS[5])
cutoff = parse(Float64, ARGS[6])
probe_row = parse(Int, ARGS[7])

isfile(campaign_file) || error("campaign file does not exist: $campaign_file")
ispath(output) && error("refusing to overwrite existing output directory: $output")
moments >= 4 || error("moments must be at least four for a compiled CUDA warm-up")
maxdim > 0 || error("maxdim must be positive")
cutoff > 0 || error("cutoff must be positive")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

# Load campaign callables before the numerical path is compiled. The legacy
# MPO construction invokes its hopping and seed closures inside TCI.
Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
case = campaign_case(Main.campaign, task)
model = case.model
model isa SquareModel || error("benchmark requires SquareModel")
nx, ny = model.size
nx == ny || error("benchmark currently requires an equal square")
N = nx * ny
levels = qtt_levels(model)
0 <= probe_row < N || error("probe row must be in 0:$(N - 1)")

representation = QTTSettings(
    encoding=:interleaved,
    tci_tol=case.representation.tci_tol,
    cutoff=cutoff,
    maxdim=maxdim,
)
parameters = legacy_parameters(model, representation, case.solver)
system = System(parameters)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)

# Match the independent vector-KPM scaling used by the small-system QTT test.
data = MPO_MeanField._kpm_data(model)
H_row = MPO_MeanField._kpm_hamiltonian(data, data.seed, zeros(length(data.bonds)))
qtt_to_row = [begin
    x, y = square_lattice_decoder(index, levels)
    x + nx * y + 1
end for index in 0:(N - 1)]
H_qtt = Matrix(H_row[qtt_to_row, qtt_to_row])
eigenvalues = eigvals(Symmetric(H_qtt))
center = (first(eigenvalues) + last(eigenvalues)) / 2
scale = (last(eigenvalues) - first(eigenvalues)) / 2 * 1.05
H_scaled = +(
    H_effective / scale,
    (-center / scale) * Identity_MPO(system.sites);
    cutoff=cutoff,
    maxdim=maxdim,
)
coefficients = MPO_MeanField._kpm_coefficients(moments, 0.0)
probe = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, probe_row)

cpu_time = @elapsed cpu_result = MPO_MeanField._qtt_mps_chebyshev_apply(
    H_scaled, probe, coefficients; cutoff=cutoff, maxdim=maxdim,
)

H_device = ITensors.adapt(CUDA.CuArray, H_scaled)
probe_device = ITensors.adapt(CUDA.CuArray, probe)
CUDA.synchronize()

# Compilation and one-time CUDA library initialization are not meaningful
# numerical timings. Warm up using a distinct short polynomial before timing.
MPO_MeanField._qtt_mps_chebyshev_apply(
    H_device, probe_device, coefficients[1:4]; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
CUDA.reclaim()
free_before = CUDA.free_memory()
gpu_time = @elapsed gpu_result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
    H_device, probe_device, coefficients; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
free_after = CUDA.free_memory()
gpu_result = ITensors.cpu(gpu_result_device.state)

cpu_amplitudes = MPO_MeanField._qtt_mps_amplitudes(cpu_result.state, system.sites)
gpu_amplitudes = MPO_MeanField._qtt_mps_amplitudes(gpu_result, system.sites)
relative_state_error = norm(cpu_amplitudes - gpu_amplitudes) / max(norm(cpu_amplitudes), sqrt(eps(Float64)))

mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()),
        "campaign" => campaign_file,
        "label" => case.label,
        "matrix_dimension" => N,
        "moments" => moments,
        "maxdim" => maxdim,
        "cutoff" => cutoff,
        "probe_row" => probe_row,
        "cpu_time_s" => cpu_time,
        "cuda_time_s" => gpu_time,
        "cuda_speedup" => cpu_time / gpu_time,
        "cpu_final_max_chi" => maxlinkdim(cpu_result.state),
        "cuda_final_max_chi" => maxlinkdim(gpu_result),
        "cpu_cuda_state_relative_error" => relative_state_error,
        "device_name" => CUDA.name(CUDA.device()),
        "device_total_memory_bytes" => CUDA.total_memory(),
        "device_free_memory_before_bytes" => free_before,
        "device_free_memory_after_bytes" => free_after,
    ))
end

println("QTT MPO--MPS CUDA benchmark complete: $output")
println("  CPU=$(cpu_time)s CUDA=$(gpu_time)s speedup=$(cpu_time / gpu_time)x")
println("  CPU/CUDA relative state error=$relative_state_error")
