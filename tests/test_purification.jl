

using CUDA
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf


include("../src/MPO_MeanField.jl")
using .MPO_MeanField

println("CUDA Functional:", CUDA.functional())
if CUDA.functional()
    println("CUDA Device Name: ", CUDA.device())
    to_gpu = cu
    to_cpu = ITensors.cpu
    cleanup = CUDA.reclaim
else
    println("CUDA is not functional. Running on CPU.")
    to_gpu = identity
    to_cpu = identity
    cleanup = () -> nothing
end

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)

# 1. System construction
# ─────────────────────────────────────────
println("\n--- Test 1: System Construction ---")


L = 20
τ = (sqrt(10)-2.0)/2.0
t(x) = -1.0 - 0.8 * cos(2π * τ * (x - 0.5))

U = 0.3
W = nothing
S(x) = 0.2 * cos(π * x)

tci_tol = 1e-6
itensors_tol = 1e-12
itensors_maxdim = 300
density = 0.5
purification_steps = 200
scf_mixing = 0.9
scf_tol = 0.1 #%
scf_max_iterations = 100

params = ModelParameters(L, t, U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys = System(params)
show(sys)

#Esta cadeira é bem melhor que a minha

#bra_cache, ket_cache = precompute_qtt_states(sys.sites)


println("\n--- Test 1a: Initial ρ₀ Trace ---")
ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
T1 = real(tr(ρ0))
Ne = round(Int, 2^L * density)
@assert abs(T1 - Ne) < 1e-6 "Trace wrong: got $T1, expected $Ne"
println("Tr(ρ₀) = $T1 ✓")



println("\n--- Test 1b: Initial ρ₀ Eigenvalues in [0,1] ---")
# All diagonal elements of the linear guess should be in [0,1]
for i in 1:1024
    val = real(MatrixChecker(ρ0, sys.sites, i, i, sys.bra_states, sys.ket_states))
    @assert -1e-6 <= val <= 1.0 + 1e-6 "Eigenvalue out of [0,1]: $val at site $i"
end
println("All diagonal elements in [0,1] ✓")



println("\n--- Test 2a: Purified ρ Idempotency ---")
ρ0 = to_gpu(ρ0)
ρ_purified = perform_purification(ρ0, params; verbose=1)
ρ_purified = to_cpu(ρ_purified) # Move back to CPU for checking

P2 = apply(ρ_purified, ρ_purified; cutoff=itensors_tol, maxdim=itensors_maxdim)
T1 = real(tr(ρ_purified))
T2 = real(tr(P2))
idem_error = abs(T1 - T2) / T1
@assert idem_error < 0.1/100 "Idempotency error too large: $(idem_error*100) %"
println("Idempotency error: $(idem_error * 100) % ✓")

println("\n--- Test 2b: Purified ρ Trace Conservation ---")
Ne_expected = round(Int, 2^L * density)
@assert abs(T1 - Ne_expected) / Ne_expected < 0.1/100 "Trace drifted: got $T1, expected $Ne_expected"
println("Tr(ρ) = $T1, Ne = $Ne_expected ✓")

println("\n--- Test 2c: Purified ρ Diagonal elements ---")
for i in 1:1024
    val = real(MatrixChecker(ρ_purified, sys.sites, i, i, sys.bra_states, sys.ket_states))
    @assert abs(val) > -1e-3 && abs(val - 1.0) > -1e-3 "Diagonal element outside [0,1]: $val at site $i"
end
println("All diagonal elements inside to 0 or 1 ✓")


# println("\n--- Test 3: Compare against Exact Diagonalization ---")
# # Build full matrix of H0
# N = 2^L
# H_matrix = zeros(ComplexF64, N, N)
# H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
# for i in 1:N
#     for j in 1:N
#         H_matrix[i,j] = MatrixChecker(H, sys.sites, i, j, sys.bra_states, sys.ket_states)
#     end
# end

# # Exact density matrix via diagonalization
# evals, evecs = eigen(Hermitian(H_matrix))
# Ne_expected = round(Int, N * density)
# ρ_exact = evecs[:, 1:Ne_expected] * evecs[:, 1:Ne_expected]'

# # Compare diagonal elements
# println("  i   ρ_exact[i,i]   ρ_purified[i,i]   error")
# max_err = 0.0
# for i in 1:N
#     val_purified = real(MatrixChecker(ρ_purified, sys.sites, i, i, sys.bra_states, sys.ket_states))
#     val_exact = real(ρ_exact[i,i])
#     err = abs(val_purified - val_exact)
#     global max_err = max(max_err, err)
#     println(@sprintf "  %2d   %8.5f        %8.5f          %8.2e" i-1 val_exact val_purified err)
# end
# @assert max_err < 1e-3 "Max diagonal error vs exact: $max_err"
# println("Max error vs exact diagonalization: $max_err ✓")