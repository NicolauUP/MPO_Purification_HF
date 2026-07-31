#!/usr/bin/env julia

"""Validate allocation-free KPM recurrences on CUDA against a CPU reference."""

using CUDA
using LinearAlgebra
using SparseArrays
using TOML

include(joinpath(@__DIR__, "kpm_local_helpers.jl"))

length(ARGS) == 1 || error(
    "usage: validate_kpm_inplace_cuda.jl OUTPUT_DIRECTORY",
)
output = abspath(only(ARGS))
ispath(output) && error("refusing to overwrite existing output: $output")
CUDA.functional() || error("CUDA is not functional on this node")
CUDA.allowscalar(false)

const N = 4096
const R = 128
const DEGREE = 256
const SEED = 510578

diagonal = collect(range(-0.55, 0.55; length=N))
offdiagonal = fill(0.1, N - 1)
hamiltonian = spdiagm(
    -1 => offdiagonal,
    0 => diagonal,
    1 => offdiagonal,
)
probes = probing_matrix(N, R, :hadamard, SEED)
coefficients = projector_coefficients(DEGREE, 0.0)

cpu_moments = kpm_trace_moments(hamiltonian, probes, DEGREE)
cpu_filtered = kpm_apply(hamiltonian, probes, coefficients)

CUDA.reclaim()
free_before = CUDA.free_memory()
gpu_hamiltonian = CUDA.CUSPARSE.CuSparseMatrixCSR(hamiltonian)
gpu_probes = CUDA.CuArray(probes)
CUDA.synchronize()
moment_calculation = @timed begin
    value = kpm_trace_moments(
        gpu_hamiltonian, gpu_probes, DEGREE;
        synchronize=CUDA.synchronize,
    )
    CUDA.synchronize()
    value
end
apply_calculation = @timed begin
    value = kpm_apply(
        gpu_hamiltonian, gpu_probes, coefficients;
        synchronize=CUDA.synchronize,
    )
    CUDA.synchronize()
    value
end
gpu_moments = moment_calculation.value
gpu_filtered = Array(apply_calculation.value)
free_after = CUDA.free_memory()

moment_error =
    norm(gpu_moments - cpu_moments) / max(norm(cpu_moments), eps())
application_error =
    norm(gpu_filtered - cpu_filtered) / max(norm(cpu_filtered), eps())
moment_error <= 1e-11 ||
    error("CUDA trace-moment relative error is $moment_error")
application_error <= 1e-11 ||
    error("CUDA KPM-application relative error is $application_error")

mkpath(output)
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "matrix_dimension" => N,
        "probes" => R,
        "degree" => DEGREE,
        "device_name" => CUDA.name(CUDA.device()),
        "moment_relative_error" => moment_error,
        "application_relative_error" => application_error,
        "moment_time_s" => moment_calculation.time,
        "application_time_s" => apply_calculation.time,
        "device_free_memory_before_bytes" => free_before,
        "device_free_memory_after_bytes" => free_after,
    ))
end
println("Allocation-free CUDA KPM validation passed: $output")
println("  moment relative error=$moment_error")
println("  application relative error=$application_error")
