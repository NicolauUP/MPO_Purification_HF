# Large-N fixed-H CSR-SP2 scaling after the N=1024 dense validation.
#
# Usage: julia --project=. scale_sparse_sp2_1d.jl OUTPUT N V2 CUTOFF
# All iterative quantities remain on the GPU. The trace and trace idempotency
# defect are reduced on device; the host receives only scalar diagnostics.
using CUDA
using LinearAlgebra
using Printf
using SparseArrays
using TOML

const TAU_AA = sqrt(2.0) - 5.0 / 6.0
const MAX_ITERATIONS = 40
const MAX_NNZ_PER_ROW = 512.0
const CONVERGENCE_TOLERANCE = 1e-8

length(ARGS) == 4 || error("usage: scale_sparse_sp2_1d.jl OUTPUT N V2 CUTOFF")
output = abspath(ARGS[1]); N = parse(Int, ARGS[2]); V2 = parse(Float64, ARGS[3]); cutoff = parse(Float64, ARGS[4])
iseven(N) || error("N must be even")
CUDA.functional() || error("CUDA is not functional on this node")
isdir(output) && error("refusing to overwrite existing output directory: $output")
mkpath(output)

function fixed_hamiltonian(N::Int, V2::Float64)
    rows = Int[]; cols = Int[]; values = Float64[]
    sizehint!(rows, 3N - 2); sizehint!(cols, 3N - 2); sizehint!(values, 3N - 2)
    for i in 1:N
        push!(rows, i); push!(cols, i); push!(values, iseven(i) ? 0.5 : -0.5)
    end
    for i in 1:N-1
        t = -1.0 - V2 * cos(2π * TAU_AA * (i + 0.5))
        push!(rows, i); push!(cols, i + 1); push!(values, t)
        push!(rows, i + 1); push!(cols, i); push!(values, t)
    end
    sparse(rows, cols, values, N, N)
end

# One GPU thread per CSR row; CSR indices from CUDA.jl are one-based here.
function csr_trace_kernel!(out, rowptr, colval, nzval, N)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if row <= N
        value = zero(eltype(nzval))
        for p in Int(rowptr[row]):(Int(rowptr[row + 1]) - 1)
            if colval[p] == row
                value = nzval[p]
                break
            end
        end
        CUDA.@atomic out[1] += value
    end
    return
end

function csr_trace(A)
    out = CUDA.zeros(Float64, 1)
    threads = 256
    @cuda threads=threads blocks=cld(size(A, 1), threads) csr_trace_kernel!(out, A.rowPtr, A.colVal, A.nzVal, size(A, 1))
    return Array(out)[1]
end

function trace_idempotency_defect(A, trace_value, particles)
    # tr(A^2)=||A||_F^2 for the real symmetric iterates in this diagnostic.
    norm_squared = sum(abs2, nonzeros(A))
    return abs(norm_squared - trace_value) / particles
end

particles = N ÷ 2
# Rigorous Gershgorin enclosure: 2 max_i |t_i| + |W_i|.
radius = 2.0 * (1.0 + abs(V2)) + 0.5
H = fixed_hamiltonian(N, V2)
P0 = sparse((radius * I - H) / (2radius))
P = CUDA.CUSPARSE.CuSparseMatrixCSR(P0)

metadata = Dict(
    "diagnostic" => "large_fixed_hamiltonian_thresholded_sparse_sp2_1d",
    "backend" => "cuda", "matrix_dimension" => N, "target_particles" => particles,
    "V2" => V2, "cutoff" => cutoff, "tau" => TAU_AA,
    "physical_staggered_potential" => "W_i = 0.5*(-1)^i", "interaction_U" => 0.0,
    "spectral_lower" => -radius, "spectral_upper" => radius,
    "maximum_iterations" => MAX_ITERATIONS, "max_nnz_per_row_guard" => MAX_NNZ_PER_ROW,
    "convergence_tolerance" => CONVERGENCE_TOLERANCE,
    "persistent_fillin_control" => "post-SpGEMM magnitude pruning",
    "preproduct_screening" => false,
    "transient_fillin_diagnostic" => "nnz_product_preprune is recorded each iteration",
    "device_name" => CUDA.name(CUDA.device()),
    "device_free_memory_before_bytes" => CUDA.free_memory(),
)
open(joinpath(output, "metadata.toml"), "w") do io; TOML.print(io, metadata); end

function run_sparse_sp2!(P, N, V2, cutoff, particles, io)
    best_score = Inf; best_iteration = 0; best_nnz = 0; best_trace = NaN; best_idempotency = NaN
    maximum_product_nnz = 0
    total_spgemm_time_s = 0.0
    termination = :maximum_iterations
    for iteration in 1:MAX_ITERATIONS
        trace_before = csr_trace(P)
        branch = trace_before >= particles ? :square : :hole
        nnz_input = nnz(P)
        CUDA.synchronize()
        spgemm_time_s = @elapsed begin product = P * P; CUDA.synchronize(); end
        nnz_product = nnz(product)
        maximum_product_nnz = max(maximum_product_nnz, nnz_product)
        total_spgemm_time_s += spgemm_time_s
        postprocess_time_s = @elapsed begin
            product = CUDA.CUSPARSE.prune(product, cutoff, 'O')
            P = branch == :square ? product : CUDA.CUSPARSE.geam(2.0, P, -1.0, product, 'O')
            P = CUDA.CUSPARSE.prune(P, cutoff, 'O')
            CUDA.synchronize()
        end
        trace_after = csr_trace(P)
        trace_error = abs(trace_after - particles) / particles
        idem = trace_idempotency_defect(P, trace_after, particles)
        output_nnz = nnz(P); per_row = output_nnz / N
        println(io, join((iteration, branch, trace_after, trace_error, idem, nnz_input, nnz_product, output_nnz, per_row, spgemm_time_s, postprocess_time_s), ',')); flush(io)
        score = max(trace_error, idem)
        if score < best_score
            best_score = score; best_iteration = iteration; best_nnz = output_nnz
            best_trace = trace_after; best_idempotency = idem
        end
        @info "sparse SP2" N V2 cutoff iteration branch trace=trace_after trace_error idem nnz=output_nnz per_row
        if per_row > MAX_NNZ_PER_ROW
            termination = :nnz_per_row_guard
            break
        elseif iteration >= 5 && score <= CONVERGENCE_TOLERANCE
            termination = :converged_trace_idempotency
            break
        elseif iteration >= 10 && score > 10 * best_score && best_score < 1e-6
            termination = :truncation_floor_rebound
            break
        end
    end
    return (; best_iteration, best_nnz, best_trace, best_idempotency,
        maximum_product_nnz, total_spgemm_time_s, termination)
end

run_summary = open(joinpath(output, "iterations.csv"), "w") do io
    println(io, "iteration,branch,trace,relative_trace_error,trace_idempotency_defect,nnz_input,nnz_product_preprune,nnz_output,nnz_per_row,spgemm_time_s,postprocess_time_s")
    run_sparse_sp2!(P, N, V2, cutoff, particles, io)
end
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, Dict(
        "termination_reason" => string(run_summary.termination), "best_iteration" => run_summary.best_iteration,
        "best_trace" => run_summary.best_trace, "best_relative_trace_error" => abs(run_summary.best_trace - particles) / particles,
        "best_trace_idempotency_defect" => run_summary.best_idempotency,
        "best_nnz" => run_summary.best_nnz, "best_nnz_per_row" => run_summary.best_nnz / N,
        "maximum_product_nnz_preprune" => run_summary.maximum_product_nnz,
        "maximum_product_nnz_per_row_preprune" => run_summary.maximum_product_nnz / N,
        "total_spgemm_time_s" => run_summary.total_spgemm_time_s,
        "device_free_memory_after_bytes" => CUDA.free_memory(),
    ))
end
