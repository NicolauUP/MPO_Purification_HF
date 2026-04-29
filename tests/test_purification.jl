# test_main.jl
__precompile__(false)
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf

include("../src/core/operators.jl")
include("../src/core/system.jl")
include("../src/hamiltonians/mpo_construction.jl")
include("../src/utils/quantics.jl")
include("src/purification/mcweeny.jl")
# include("src/tci/modulations.jl")

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)

# 1. System construction
# ─────────────────────────────────────────
println("\n--- Test 1: System Construction ---")
L = 6
t = -1.0
U = 0.0
W(x) = 0.5 * cos(π * x)
S = nothing
tci_tol = 1e-6
itensors_tol = 1e-10
itensors_maxdim = 100
density = 0.5
purification_steps = 40
scf_mixing = 0.5
scf_tol = 0.1 #%
scf_max_iterations = 100

params = ModelParameters(L, t, U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys = System(params)
println("System constructed successfully:")
println()
show(sys)

#Esta cadeira é bem melhor que a minha

bra_cache, ket_cache = precompute_qtt_states(sys.sites)


println("\n--- Test 1a: Initial ρ₀ Trace ---")
ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
T1 = real(tr(ρ0))
Ne_expected = round(Int, 2^L * density)
@assert abs(T1 - Ne_expected) < 1e-6 "Trace wrong: got $T1, expected $Ne_expected"
println("Tr(ρ₀) = $T1 ✓")


println("\n--- Test 1b: Initial ρ₀ Eigenvalues in [0,1] ---")
# All diagonal elements of the linear guess should be in [0,1]
for i in 1:2^L
    val = real(MatrixChecker(ρ0, sys.sites, i, i, bra_cache, ket_cache))
    @assert -1e-6 <= val <= 1.0 + 1e-6 "Eigenvalue out of [0,1]: $val at site $i"
end
println("All diagonal elements in [0,1] ✓")



println("\n--- Test 2a: Purified ρ Idempotency ---")
ρ_purified = perform_purification(ρ0, params; verbose=1)


P2 = apply(ρ_purified, ρ_purified; cutoff=itensors_tol, maxdim=itensors_maxdim)
T1 = real(tr(ρ_purified))
T2 = real(tr(P2))
idem_error = abs(T1 - T2) / T1
@assert idem_error < 0.1/100 "Idempotency error too large: $(idem_error*100) %"
println("Idempotency error: $(idem_error * 100) % ✓")

println("\n--- Test 2b: Purified ρ Trace Conservation ---")
Ne_expected = round(Int, 2^L * density)
@assert abs(T1 - Ne_expected) < 0.1/100 "Trace drifted: got $T1, expected $Ne_expected"
println("Tr(ρ) = $T1, Ne = $Ne_expected ✓")

println("\n--- Test 2c: Purified ρ Eigenvalues in {0,1} ---")
# Diagonal elements should be close to 0 or 1
for i in 1:2^L
    val = real(MatrixChecker(ρ_purified, sys.sites, i, i, bra_cache, ket_cache))
    @assert val < 1e-3 || val > 1.0 - 1e-3 "Eigenvalue not close to 0 or 1: $val at $i"
end
println("All diagonal elements close to 0 or 1 ✓")


println("\n--- Test 3: Compare against Exact Diagonalization ---")
# Build full matrix of H0
N = 2^L
H_matrix = zeros(ComplexF64, N, N)
for i in 1:N
    for j in 1:N
        H_matrix[i,j] = MatrixChecker(H0, sys.sites, i, j, bra_cache, ket_cache)
    end
end

# Exact density matrix via diagonalization
evals, evecs = eigen(Hermitian(H_matrix))
Ne_expected = round(Int, N * density)
ρ_exact = evecs[:, 1:Ne_expected] * evecs[:, 1:Ne_expected]'

# Compare diagonal elements
println("  i   ρ_exact[i,i]   ρ_purified[i,i]   error")
max_err = 0.0
for i in 1:N
    val_purified = real(MatrixChecker(ρ_purified, sys.sites, i, i, bra_cache, ket_cache))
    val_exact = real(ρ_exact[i,i])
    err = abs(val_purified - val_exact)
    global max_err = max(max_err, err)
    println(@sprintf "  %2d   %8.5f        %8.5f          %8.2e" i-1 val_exact val_purified err)
end
@assert max_err < 1e-3 "Max diagonal error vs exact: $max_err"
println("Max error vs exact diagonalization: $max_err ✓")