"""Validate thresholded CSR-SP2 against an exact 1D projector.

This is deliberately independent of the MPO and SCF implementations.  It
uses a real-space open chain with quasiperiodic hopping and a *physical*
staggered onsite field.  The latter opens a gap while U=0 keeps the experiment
strictly fixed-H.

Usage
-----
    julia --project=. benchmark/sparse_sp2_gpu/diagnose_sparse_sp2_1d.jl OUTPUT [cuda|cpu] [N]

The output contains one row per (V2, cutoff, SP2 iteration).  At N=1024 an
exact dense projector is retained exclusively as a validation reference.
"""

using CUDA
using LinearAlgebra
using Printf
using SparseArrays
using TOML

const TAU_AA = sqrt(2.0) - 5.0 / 6.0
const V2_VALUES = (0.0, 0.5, 2.0)
const NEAR_MACHINE_CUTOFF = 64 * eps(Float64)
const CUTOFFS = (NEAR_MACHINE_CUTOFF, 1e-12, 1e-10, 1e-8, 1e-6)
const MAXIMUM_ITERATIONS = 32

length(ARGS) >= 1 || error("usage: diagnose_sparse_sp2_1d.jl OUTPUT [cuda|cpu] [N]")
output = abspath(ARGS[1])
backend = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :cuda
N = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1024
N > 1 || error("N must exceed one")
iseven(N) || error("N must be even at half filling")
backend in (:cuda, :cpu) || error("backend must be cuda or cpu")

isdir(output) && error("refusing to overwrite existing output directory: $output")
mkpath(output)

"Open-chain hopping Hamiltonian with W_i = 0.5*(-1)^i."
function fixed_hamiltonian(N::Int, V2::Float64)
    rows = Int[]
    cols = Int[]
    values = Float64[]
    sizehint!(rows, 3N - 2)
    sizehint!(cols, 3N - 2)
    sizehint!(values, 3N - 2)
    for i in 1:N
        push!(rows, i); push!(cols, i); push!(values, iseven(i) ? 0.5 : -0.5)
    end
    for i in 1:(N - 1)
        hopping = -1.0 - V2 * cos(2π * TAU_AA * (i + 0.5))
        push!(rows, i); push!(cols, i + 1); push!(values, hopping)
        push!(rows, i + 1); push!(cols, i); push!(values, hopping)
    end
    return sparse(rows, cols, values, N, N)
end

function exact_projector(H::SparseMatrixCSC{Float64,Int}, particles::Int)
    spectral = eigen(Symmetric(Matrix(H)))
    vectors = spectral.vectors[:, 1:particles]
    return vectors * vectors', spectral.values
end

function initial_sp2_matrix(H::SparseMatrixCSC{Float64,Int}, lower::Float64, upper::Float64)
    # Low energies map to one and high energies to zero.
    return sparse((upper * I - H) / (upper - lower))
end

function sparse_prune_cpu(A::SparseMatrixCSC{Float64,Int}, cutoff::Float64)
    rows, cols, values = findnz(A)
    keep = abs.(values) .> cutoff
    return sparse(rows[keep], cols[keep], values[keep], size(A)...)
end

function sparse_symmetrize_cpu(A::SparseMatrixCSC{Float64,Int}, cutoff::Float64)
    return sparse_prune_cpu(sparse(0.5 * (A + A')), cutoff)
end

function sparse_prune_gpu(A, cutoff::Float64)
    # CUDA.jl binds cuSPARSE's CSR->CSR magnitude pruning routine.  Index 'O'
    # denotes one-based CSR indices, matching CUDA.jl's constructors.
    return CUDA.CUSPARSE.prune(A, cutoff, 'O')
end

function gpu_spgemm(A, B)
    # This dispatches to generic cuSPARSE SpGEMM.  The output is CSR and sorted.
    return A * B
end

function snapshot_metrics(P_sparse::SparseMatrixCSC{Float64,Int}, P_exact::Matrix{Float64})
    # Avoid a dense O(N^3) validation product.  The exact reference is dense at
    # N=1024, but all diagnostics of the sparse iterate are evaluated from its
    # stored entries.  Since P_exact is a rank-Ne orthogonal projector,
    # ||P-P_exact||_F^2 = ||P||_F^2 + Ne - 2<P,P_exact>.
    N = size(P_sparse, 1)
    particles = N ÷ 2
    rows, cols, values = findnz(P_sparse)
    norm_squared = sum(abs2, values)
    overlap = sum(values[k] * P_exact[rows[k], cols[k]] for k in eachindex(values))
    relative_error_squared = max(0.0, (norm_squared + particles - 2.0 * overlap) / particles)
    rho_error = [P_sparse[i, i] - P_exact[i, i] for i in 1:N]
    bond_error = [P_sparse[i, i + 1] - P_exact[i, i + 1] for i in 1:(N - 1)]
    return (
        trace=sum(P_sparse[i, i] for i in 1:N),
        relative_projector_error=sqrt(relative_error_squared),
        hermiticity_residual=norm(P_sparse - P_sparse') / max(sqrt(norm_squared), eps()),
        # For a Hermitian P, tr(P^2) = ||P||_F^2. This is the inexpensive
        # trace-normalized idempotency diagnostic used by the main solver.
        trace_idempotency_defect=abs(norm_squared - sum(P_sparse[i, i] for i in 1:N)) / particles,
        density_max_abs_error=maximum(abs, rho_error),
        density_rms_error=sqrt(sum(abs2, rho_error) / N),
        bond_max_abs_error=maximum(abs, bond_error),
        bond_rms_error=sqrt(sum(abs2, bond_error) / length(bond_error)),
    )
end

function write_header(io)
    println(io, join((
        "v2", "cutoff", "iteration", "branch", "nnz_input", "nnz_product_preprune",
        "nnz_output", "nnz_per_row", "spgemm_time_s", "postprocess_time_s",
        "trace", "relative_projector_error", "hermiticity_residual",
        "trace_idempotency_defect",
        "density_max_abs_error", "density_rms_error", "bond_max_abs_error", "bond_rms_error",
    ), ','))
end

function write_row(io; v2, cutoff, iteration, branch, nnz_input, nnz_product_preprune,
                   nnz_output, nnz_per_row, spgemm_time_s, postprocess_time_s, metrics)
    println(io, join((
        v2, cutoff, iteration, branch, nnz_input, nnz_product_preprune,
        nnz_output, nnz_per_row, spgemm_time_s, postprocess_time_s, metrics.trace,
        metrics.relative_projector_error, metrics.hermiticity_residual,
        metrics.trace_idempotency_defect,
        metrics.density_max_abs_error, metrics.density_rms_error,
        metrics.bond_max_abs_error, metrics.bond_rms_error,
    ), ','))
end

particles = N ÷ 2
metadata = Dict(
    "diagnostic" => "fixed_hamiltonian_thresholded_sparse_sp2_1d",
    "backend" => string(backend),
    "matrix_dimension" => N,
    "target_particles" => particles,
    "interaction_U" => 0.0,
    "physical_staggered_potential" => "W_i = 0.5*(-1)^i",
    "hopping" => "t_i = -1 - V2*cos(2*pi*tau*(i+1/2))",
    "tau" => TAU_AA,
    "v2_values" => collect(V2_VALUES),
    "cutoffs" => collect(CUTOFFS),
    "machine_precision_cutoff" => NEAR_MACHINE_CUTOFF,
    "sp2_idempotency_tolerance" => 2e-8,
    "sp2_relative_trace_tolerance" => 2e-8,
    "maximum_iterations" => MAXIMUM_ITERATIONS,
    "cuda_functional" => CUDA.functional(),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

if backend == :cuda
    CUDA.functional() || error("CUDA is not functional on this node")
end

open(joinpath(output, "iterations.csv"), "w") do io
    write_header(io)
    for V2 in V2_VALUES
        @info "fixed-H sparse SP2" N V2 backend
        H = fixed_hamiltonian(N, V2)
        P_exact, eigenvalues = exact_projector(H, particles)
        lower, upper = extrema(eigenvalues)
        initial = initial_sp2_matrix(H, lower, upper)

        for cutoff in CUTOFFS
            P_cpu = copy(initial)
            P_gpu = backend == :cuda ? CUDA.CUSPARSE.CuSparseMatrixCSR(P_cpu) : nothing
            for iteration in 1:MAXIMUM_ITERATIONS
                nnz_input = nnz(P_cpu)
                product, spgemm_time_s = if backend == :cuda
                    CUDA.synchronize()
                    elapsed = @elapsed begin
                        product_gpu = gpu_spgemm(P_gpu, P_gpu)
                        CUDA.synchronize()
                    end
                    product_gpu, elapsed
                else
                    elapsed = @elapsed product_cpu = P_cpu * P_cpu
                    product_cpu, elapsed
                end
                nnz_product = nnz(product)

                branch = tr(P_cpu) >= particles ? :square : :hole
                next_cpu, postprocess_time_s = if backend == :cuda
                    elapsed = @elapsed begin
                        product_pruned_gpu = sparse_prune_gpu(product, cutoff)
                        candidate_gpu = branch == :square ? product_pruned_gpu :
                            CUDA.CUSPARSE.geam(2.0, P_gpu, -1.0, product_pruned_gpu, 'O')
                        candidate_gpu = sparse_prune_gpu(candidate_gpu, cutoff)
                        CUDA.synchronize()
                    end
                    # A host copy is made only for validation and branch selection.
                    # The SP2 multiply, affine update, and magnitude pruning all ran
                    # on the H100 above.
                    P_gpu = candidate_gpu
                    SparseMatrixCSC(Array(P_gpu)), elapsed
                else
                    elapsed = @elapsed begin
                        candidate = branch == :square ? product : sparse(2.0 * P_cpu - product)
                        candidate = sparse_symmetrize_cpu(candidate, cutoff)
                    end
                    candidate, elapsed
                end
                P_cpu = next_cpu
                metrics = snapshot_metrics(P_cpu, P_exact)
                write_row(io;
                    v2=V2, cutoff=cutoff, iteration=iteration, branch=branch,
                    nnz_input=nnz_input, nnz_product_preprune=nnz_product,
                    nnz_output=nnz(P_cpu), nnz_per_row=nnz(P_cpu) / N,
                    spgemm_time_s=spgemm_time_s, postprocess_time_s=postprocess_time_s,
                    metrics=metrics,
                )
                flush(io)
                @info "SP2" V2 cutoff iteration branch trace=metrics.trace nnz=nnz(P_cpu) density_rms_error=metrics.density_rms_error
                metrics.relative_projector_error < 2e-8 &&
                    abs(metrics.trace - particles) / particles < 2e-8 && break
            end
        end
    end
end
