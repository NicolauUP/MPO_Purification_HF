#!/usr/bin/env julia

"""CPU/CUDA QTT MPO--MPS KPM comparison for a legacy fixed-H campaign.

This driver deliberately avoids dense matrices and does exactly one rank-one
Hadamard probe. It is the scale-relevant backend test before using CUDA for a
multi-probe density oracle. The campaign must supply a symmetric fixed-H
spectral interval and must be particle-hole symmetric at μ=0.
"""

using Dates
using LinearAlgebra
using TOML
using MPO_MeanField
using ITensors, ITensorMPS
using CUDA

length(ARGS) == 7 || error(
    "usage: $(PROGRAM_FILE) LEGACY_CAMPAIGN.jl TASK OUTPUT MOMENTS MAXDIM CUTOFF PROBE_ROW",
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
moments >= 4 || error("moments must be at least four")
maxdim > 0 || error("maxdim must be positive")
cutoff > 0 || error("cutoff must be positive")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

Base.include(Main, campaign_file)
isdefined(Main, :campaign) || error("campaign file did not define `campaign`")
spec = Main.campaign.runs[task]
params = spec.params
params isa ParametersSquare || error("benchmark requires ParametersSquare")
N = 1 << params.L
0 <= probe_row < N || error("probe row must be in 0:$(N - 1)")
lower, upper = validate_spectral_bounds(spec.spectral_bounds...)
abs(lower + upper) <= 1e-12 || error("large QTT CUDA benchmark requires symmetric spectral bounds")

effective_params = ParametersSquare(
    L=params.L, t=params.t, U=params.U, W=params.W, S=params.S,
    tci_tol=params.tci_tol, itensors_tol=cutoff, itensors_maxdim=maxdim,
    density=params.density, purification_steps=params.purification_steps,
    scf_mixing=params.scf_mixing, scf_tol=params.scf_tol,
    scf_max_iterations=params.scf_max_iterations,
)
system = System(effective_params)
H_effective = +(system.H0, system.VH, system.VF; cutoff=cutoff, maxdim=maxdim)
H_scaled = H_effective / ((upper - lower) / 2)
coefficients = MPO_MeanField._kpm_coefficients(moments, 0.0)
probe = MPO_MeanField._qtt_hadamard_probe_mps(system.sites, probe_row)

cpu_time = @elapsed cpu_result = MPO_MeanField._qtt_mps_chebyshev_apply(
    H_scaled, probe, coefficients; cutoff=cutoff, maxdim=maxdim,
)

H_device = ITensors.adapt(CUDA.CuArray, H_scaled)
probe_device = ITensors.adapt(CUDA.CuArray, probe)
CUDA.synchronize()
MPO_MeanField._qtt_mps_chebyshev_apply(
    H_device, probe_device, coefficients[1:4]; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
CUDA.reclaim()
free_before = CUDA.free_memory()
cuda_time = @elapsed cuda_result_device = MPO_MeanField._qtt_mps_chebyshev_apply(
    H_device, probe_device, coefficients; cutoff=cutoff, maxdim=maxdim,
)
CUDA.synchronize()
free_after = CUDA.free_memory()
cuda_result = ITensors.cpu(cuda_result_device.state)

# Compare complete MPS states without materialising N=65536 amplitudes.
difference = +(cpu_result.state, -cuda_result; cutoff=cutoff, maxdim=maxdim)
relative_state_error = norm(difference) / max(norm(cpu_result.state), sqrt(eps(Float64)))

mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "created_at" => string(now()), "campaign" => campaign_file,
        "label" => spec.label, "matrix_dimension" => N, "moments" => moments,
        "maxdim" => maxdim, "cutoff" => cutoff, "probe_row" => probe_row,
        "spectral_lower" => lower, "spectral_upper" => upper,
        "chemical_potential" => 0.0, "cpu_time_s" => cpu_time,
        "cuda_time_s" => cuda_time, "cuda_speedup" => cpu_time / cuda_time,
        "cpu_final_max_chi" => maxlinkdim(cpu_result.state),
        "cuda_final_max_chi" => maxlinkdim(cuda_result),
        "cpu_cuda_state_relative_error" => relative_state_error,
        "device_name" => CUDA.name(CUDA.device()),
        "device_total_memory_bytes" => CUDA.total_memory(),
        "device_free_memory_before_bytes" => free_before,
        "device_free_memory_after_bytes" => free_after,
    ))
end

println("large QTT MPO--MPS CUDA benchmark complete: $output")
println("  CPU=$(cpu_time)s CUDA=$(cuda_time)s speedup=$(cpu_time / cuda_time)x")
println("  CPU/CUDA relative state error=$relative_state_error")
