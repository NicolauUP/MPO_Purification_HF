"""Thresholded CSR-SP2 Hartree--Fock SCF for the open 1D hopping Aubry--Andre chain.

Usage:
  julia --project=. run_sparse_sp2_scf_1d.jl OUTPUT N V2 U SEED CUTOFF [--validate-dense]

The physical convention intentionally matches the existing 1D MPO path:
  H0[i,i+1] = -1 - V2*cos(2*pi*tau*(i+1/2)),
  VH[i,i]    = U*(n[i-1] + n[i+1]),
  VF[i,i+1]  = -U*real(rho[i,i+1]).

The seed is an initial staggered Hartree field, not an external potential.
All SP2 products, pruning, trace evaluation, and local-density extraction run
on the GPU. The O(N) mean-field update is deliberately host-side in this first
diagnostic, making the sparse purification cost separable from field handling.
"""

using CUDA
using LinearAlgebra
using SparseArrays
using TOML

const TAU_AA = sqrt(2.0) - 5.0 / 6.0
const MAX_SCF_ITERATIONS = parse(Int, get(ENV, "SPARSE_MAX_SCF_ITERATIONS", "30"))
const MAX_SP2_ITERATIONS = 60
const SCF_TOLERANCE = 1e-4
const SP2_RELATIVE_TRACE_TOLERANCE = 1e-6
const SP2_IDEMPOTENCY_TOLERANCE = 2e-4
const STABLE_ITERATIONS = 3
const MIXING = 0.5
const MAX_NNZ_PER_ROW = 512.0

length(ARGS) in (6, 7) || error(
    "usage: run_sparse_sp2_scf_1d.jl OUTPUT N V2 U SEED CUTOFF [--validate-dense]",
)
output = abspath(ARGS[1])
N = parse(Int, ARGS[2])
V2 = parse(Float64, ARGS[3])
U = parse(Float64, ARGS[4])
seed_amplitude = parse(Float64, ARGS[5])
cutoff = parse(Float64, ARGS[6])
validate_dense = length(ARGS) == 7 && ARGS[7] == "--validate-dense"
length(ARGS) == 7 && !validate_dense && error("unrecognized option $(ARGS[7])")

iseven(N) || error("N must be even at half filling")
N > 2 || error("N must exceed two")
U >= 0 || error("this diagnostic currently assumes repulsive U >= 0")
CUDA.functional() || error("CUDA is not functional on this node")
isdir(output) && error("refusing to overwrite existing output directory: $output")
validate_dense && N > 2048 && error("dense validation is restricted to N <= 2048")
mkpath(output)

"Build H0 + VH + VF as an open-chain sparse matrix on the host."
function sparse_effective_hamiltonian(
    N::Int, V2::Float64, hartree::Vector{Float64}, fock::Vector{Float64},
)
    length(hartree) == N || error("Hartree field length mismatch")
    length(fock) == N - 1 || error("Fock field length mismatch")
    rows = Int[]; columns = Int[]; values = Float64[]
    sizehint!(rows, 3N - 2); sizehint!(columns, 3N - 2); sizehint!(values, 3N - 2)
    for i in 1:N
        push!(rows, i); push!(columns, i); push!(values, hartree[i])
    end
    for i in 1:(N - 1)
        hopping = -1.0 - V2 * cos(2π * TAU_AA * (i + 0.5)) + fock[i]
        push!(rows, i); push!(columns, i + 1); push!(values, hopping)
        push!(rows, i + 1); push!(columns, i); push!(values, hopping)
    end
    return sparse(rows, columns, values, N, N)
end

"A symmetric Gershgorin enclosure for the real tridiagonal effective Hamiltonian."
function spectral_radius_bound(hartree::Vector{Float64}, hamiltonian::SparseMatrixCSC)
    N = length(hartree)
    bound = 0.0
    for i in 1:N
        row_sum = abs(hartree[i])
        i > 1 && (row_sum += abs(hamiltonian[i, i - 1]))
        i < N && (row_sum += abs(hamiltonian[i, i + 1]))
        bound = max(bound, row_sum)
    end
    return max(bound, 1.0)
end

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
    @cuda threads=threads blocks=cld(size(A, 1), threads) csr_trace_kernel!(
        out, A.rowPtr, A.colVal, A.nzVal, size(A, 1),
    )
    return Array(out)[1]
end

"Extract only n_i and real rho_{i,i+1}; no sparse matrix is copied to host."
function csr_local_observables_kernel!(density, bonds, rowptr, colval, nzval, N)
    row = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if row <= N
        diagonal = zero(eltype(nzval))
        right_bond = zero(eltype(nzval))
        for p in Int(rowptr[row]):(Int(rowptr[row + 1]) - 1)
            column = colval[p]
            column == row && (diagonal = nzval[p])
            row < N && column == row + 1 && (right_bond = nzval[p])
        end
        density[row] = diagonal
        row < N && (bonds[row] = right_bond)
    end
    return
end

function local_observables(A, N::Int)
    density = CUDA.zeros(Float64, N)
    bonds = CUDA.zeros(Float64, N - 1)
    threads = 256
    @cuda threads=threads blocks=cld(N, threads) csr_local_observables_kernel!(
        density, bonds, A.rowPtr, A.colVal, A.nzVal, N,
    )
    return Array(density), Array(bonds)
end

function trace_idempotency_defect(A, trace_value, particles)
    # The sparse iterates are real symmetric up to cuSPARSE roundoff, so
    # tr(P^2)=||P||_F^2 is the inexpensive idempotency diagnostic.
    return abs(sum(abs2, nonzeros(A)) - trace_value) / particles
end

function sparse_sp2(
    hamiltonian::SparseMatrixCSC{Float64,Int}, radius::Float64, particles::Int,
    cutoff::Float64, scf_iteration::Int, history_io::IO;
    preproduct_cutoff::Float64=0.0,
)
    initial = sparse((radius * I - hamiltonian) / (2radius))
    P = CUDA.CUSPARSE.CuSparseMatrixCSR(initial)
    best_score = Inf
    best = nothing
    termination = :maximum_iterations
    max_product_nnz = 0

    for iteration in 1:MAX_SP2_ITERATIONS
        trace_before = csr_trace(P)
        branch = trace_before >= particles ? :square : :hole
        CUDA.synchronize()
        preprune_time = @elapsed begin
            if preproduct_cutoff > cutoff
                P = CUDA.CUSPARSE.prune(P, preproduct_cutoff, 'O')
            end
            CUDA.synchronize()
        end
        nnz_input = nnz(P)
        CUDA.synchronize()
        product_time = @elapsed begin
            product = P * P
            CUDA.synchronize()
        end
        nnz_product = nnz(product)
        max_product_nnz = max(max_product_nnz, nnz_product)
        postprocess_time = @elapsed begin
            product = CUDA.CUSPARSE.prune(product, cutoff, 'O')
            P = branch == :square ? product : CUDA.CUSPARSE.geam(2.0, P, -1.0, product, 'O')
            P = CUDA.CUSPARSE.prune(P, cutoff, 'O')
            CUDA.synchronize()
        end
        trace_after = csr_trace(P)
        trace_error = abs(trace_after - particles) / particles
        idempotency = trace_idempotency_defect(P, trace_after, particles)
        output_nnz = nnz(P)
        per_row = output_nnz / size(P, 1)
        println(history_io, join((
            scf_iteration, iteration, branch, trace_after, trace_error, idempotency,
            nnz_input, nnz_product, output_nnz, per_row, preprune_time,
            product_time, postprocess_time,
        ), ','))
        flush(history_io)
        score = max(trace_error, idempotency)
        if score < best_score
            best_score = score
            best = (trace=trace_after, trace_error, idempotency, iteration, output_nnz, per_row)
        end
        if per_row > MAX_NNZ_PER_ROW
            termination = :nnz_per_row_guard
            break
        elseif trace_error <= SP2_RELATIVE_TRACE_TOLERANCE && idempotency <= SP2_IDEMPOTENCY_TOLERANCE
            termination = :converged_trace_idempotency
            break
        end
    end
    return P, merge(best, (; termination, max_product_nnz))
end

function mean_field_candidates(density::Vector{Float64}, bonds::Vector{Float64}, U::Float64)
    N = length(density)
    hartree = zeros(Float64, N)
    for i in 1:N
        i > 1 && (hartree[i] += U * density[i - 1])
        i < N && (hartree[i] += U * density[i + 1])
    end
    fock = -U .* bonds
    return hartree, fock
end

function field_residual(candidate_hartree, candidate_fock, hartree, fock)
    return max(
        maximum(abs, candidate_hartree .- hartree),
        maximum(abs, candidate_fock .- fock),
    )
end

function dense_scf_reference(N, V2, U, seed_amplitude)
    particles = N ÷ 2
    hartree = [iseven(i) ? seed_amplitude : -seed_amplitude for i in 1:N]
    fock = zeros(Float64, N - 1)
    stable = 0
    for iteration in 1:MAX_SCF_ITERATIONS
        H = Matrix(sparse_effective_hamiltonian(N, V2, hartree, fock))
        eigensystem = eigen(Symmetric(H))
        occupied = eigensystem.vectors[:, 1:particles]
        rho = occupied * occupied'
        density = real.(diag(rho))
        bonds = [real(rho[i, i + 1]) for i in 1:(N - 1)]
        candidate_hartree, candidate_fock = mean_field_candidates(density, bonds, U)
        residual = field_residual(candidate_hartree, candidate_fock, hartree, fock)
        stable = residual <= SCF_TOLERANCE ? stable + 1 : 0
        stable >= STABLE_ITERATIONS && return (; converged=true, iteration, density, bonds)
        hartree .= MIXING .* candidate_hartree .+ (1 - MIXING) .* hartree
        fock .= MIXING .* candidate_fock .+ (1 - MIXING) .* fock
    end
    return (; converged=false, iteration=MAX_SCF_ITERATIONS, density=zeros(N), bonds=zeros(N - 1))
end

particles = N ÷ 2
metadata = Dict(
    "diagnostic" => "thresholded_sparse_sp2_hartree_fock_scf_1d",
    "backend" => "cuda", "matrix_dimension" => N, "target_particles" => particles,
    "V2" => V2, "interaction_U" => U, "seed_amplitude" => seed_amplitude,
    "seed_definition" => "initial Hartree field S_i = seed_amplitude*(-1)^i with + sign on even i",
    "hopping" => "t_i = -1 - V2*cos(2*pi*tau*(i+1/2))",
    "hartree_convention" => "VH_i = U*(n_(i-1)+n_(i+1))",
    "fock_convention" => "VF_(i,i+1) = -U*real(rho_(i,i+1))",
    "tau" => TAU_AA, "cutoff" => cutoff,
    "scf_tolerance" => SCF_TOLERANCE, "scf_mixing" => MIXING,
    "max_scf_iterations" => MAX_SCF_ITERATIONS,
    "sp2_relative_trace_tolerance" => SP2_RELATIVE_TRACE_TOLERANCE,
    "sp2_idempotency_tolerance" => SP2_IDEMPOTENCY_TOLERANCE,
    "max_sp2_iterations" => MAX_SP2_ITERATIONS,
    "max_nnz_per_row_guard" => MAX_NNZ_PER_ROW,
    "preproduct_screening" => haskey(ENV, "SPARSE_PREPRODUCT_CUTOFF"),
    "preproduct_cutoff" => parse(Float64, get(ENV, "SPARSE_PREPRODUCT_CUTOFF", "0.0")),
    "device_name" => CUDA.name(CUDA.device()),
    "device_free_memory_before_bytes" => CUDA.free_memory(),
)
open(joinpath(output, "metadata.toml"), "w") do io
    TOML.print(io, metadata)
end

let
hartree = [iseven(i) ? seed_amplitude : -seed_amplitude for i in 1:N]
fock = zeros(Float64, N - 1)
previous_density = nothing
stable = 0
converged = false
termination = :maximum_scf_iterations
final_density = zeros(Float64, N)
final_bonds = zeros(Float64, N - 1)
final_sp2 = nothing

open(joinpath(output, "sp2_history.csv"), "w") do sp2_io
    println(sp2_io, "scf_iteration,sp2_iteration,branch,trace,relative_trace_error,trace_idempotency_defect,nnz_input,nnz_product_preprune,nnz_output,nnz_per_row,operand_prune_time_s,spgemm_time_s,postprocess_time_s")
    open(joinpath(output, "scf_history.csv"), "w") do scf_io
        println(scf_io, "iteration,spectral_radius,trace,relative_trace_error,trace_idempotency_defect,sp2_iterations,sp2_termination,sp2_nnz_per_row,sp2_product_nnz_per_row,field_residual,density_change_max_abs,stable_iterations")
        for scf_iteration in 1:MAX_SCF_ITERATIONS
            H = sparse_effective_hamiltonian(N, V2, hartree, fock)
            radius = spectral_radius_bound(hartree, H)
            P, sp2 = sparse_sp2(H, radius, particles, cutoff, scf_iteration, sp2_io;
                preproduct_cutoff=parse(Float64, get(ENV, "SPARSE_PREPRODUCT_CUTOFF", "0.0")))
            sp2.termination == :nnz_per_row_guard && begin
                termination = :sparse_nnz_guard
                final_sp2 = sp2
                break
            end
            density, bonds = local_observables(P, N)
            candidate_hartree, candidate_fock = mean_field_candidates(density, bonds, U)
            residual = field_residual(candidate_hartree, candidate_fock, hartree, fock)
            density_change = isnothing(previous_density) ? Inf : maximum(abs, density .- previous_density)
            stable = residual <= SCF_TOLERANCE ? stable + 1 : 0
            println(scf_io, join((
                scf_iteration, radius, sp2.trace, sp2.trace_error, sp2.idempotency,
                sp2.iteration, sp2.termination, sp2.per_row, sp2.max_product_nnz / N,
                residual, density_change, stable,
            ), ','))
            flush(scf_io)
            @info "sparse SP2 SCF iteration=$scf_iteration trace=$(sp2.trace) idem=$(sp2.idempotency) sp2_nnz_per_row=$(sp2.per_row) field_residual=$residual density_change=$density_change"
            final_density .= density
            final_bonds .= bonds
            final_sp2 = sp2
            if stable >= STABLE_ITERATIONS
                converged = true
                termination = :converged
                break
            end
            hartree .= MIXING .* candidate_hartree .+ (1 - MIXING) .* hartree
            fock .= MIXING .* candidate_fock .+ (1 - MIXING) .* fock
            previous_density = density
        end
    end
end

open(joinpath(output, "observables.toml"), "w") do io
    TOML.print(io, Dict(
        "scf_converged" => converged,
        "scf_termination_reason" => string(termination),
        "particle_number" => isnothing(final_sp2) ? NaN : final_sp2.trace,
        "trace_error" => isnothing(final_sp2) ? NaN : final_sp2.trace_error,
        "trace_idempotency_defect" => isnothing(final_sp2) ? NaN : final_sp2.idempotency,
        "final_nnz_per_row" => isnothing(final_sp2) ? NaN : final_sp2.per_row,
        "peak_product_nnz_per_row" => isnothing(final_sp2) ? NaN : final_sp2.max_product_nnz / N,
        "device_free_memory_after_bytes" => CUDA.free_memory(),
    ))
end

open(joinpath(output, "density.csv"), "w") do io
    println(io, "site,density")
    for i in 1:N
        println(io, "$(i),$(final_density[i])")
    end
end
open(joinpath(output, "bond_order.csv"), "w") do io
    println(io, "left_site,right_site,bond_order")
    for i in 1:(N - 1)
        println(io, "$(i),$(i + 1),$(final_bonds[i])")
    end
end

if validate_dense
    dense = dense_scf_reference(N, V2, U, seed_amplitude)
    open(joinpath(output, "dense_validation.toml"), "w") do io
        TOML.print(io, Dict(
            "dense_scf_converged" => dense.converged,
            "dense_scf_iterations" => dense.iteration,
            "density_max_abs_error" => maximum(abs, final_density .- dense.density),
            "density_rms_error" => sqrt(sum(abs2, final_density .- dense.density) / N),
            "bond_max_abs_error" => maximum(abs, final_bonds .- dense.bonds),
            "bond_rms_error" => sqrt(sum(abs2, final_bonds .- dense.bonds) / (N - 1)),
        ))
    end
end

println("Sparse SP2 SCF complete: $output")
println("  converged=$converged termination=$termination")
end
