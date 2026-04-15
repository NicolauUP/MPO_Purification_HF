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

# 2. MPO check
println("\n--- Test 2: MPO Check ---")
H0 = sys.H0
for i in 1:2^L
    for j in 1:2^L
        val = MatrixChecker(H0, sys.sites, i , j , bra_cache, ket_cache)
        if abs(val) > 1e-6
            print(@sprintf "<%.0f|H0|%.0f> = %8.3f   " i-1 j-1 val)
        end
    end
    println()
end



# 3. Purification Check

ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
println("\n--- Test 3: Purification Check ---")
ρ_purified = perform_purification(ρ0, params; verbose=1)


for i in 1:2^L
    @time val = MatrixChecker(ρ_purified, sys.sites, i , i ,bra_cache, ket_cache) # Check diagonal elements (expecting values close to 0 or 1)

    println(@sprintf "<%.0f|ρ|%.0f> = %8.3f   " i - 1 i - 1 val)
end


for i in 1:2^L
    @time val = MatrixChecker(ρ_purified, sys.sites, i , i) # Check diagonal elements (expecting values close to 0 or 1)

    println(@sprintf "<%.0f|ρ|%.0f> = %8.3f   " i - 1 i - 1 val)
end
