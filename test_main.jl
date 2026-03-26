# test_main.jl
__precompile__(false)
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf

include("src/core/operators.jl")
include("src/core/system.jl")
include("src/hamiltonians/mpo_construction.jl")
include("src/utils/quantics.jl")
# include("src/purification/mcweeny.jl")
# include("src/tci/modulations.jl")

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)

# 1. System construction
# ─────────────────────────────────────────
println("\n--- Test 1: System Construction ---")
L = 4
t(x) = 1.0 * cos(pi * √5 * (x+0.5))
U = 0.0
W(x) = 0.5 * cos(π * x) 
tci_tol = 1e-6
itensors_tol = 1e-10
itensors_maxdim = 100
params = ModelParameters(L, t, U, W, tci_tol, itensors_tol, itensors_maxdim)
sys = System(params)
println("System constructed successfully:")
println()
show(sys)



# 2. MPO check
println("\n--- Test 2: MPO Check ---")
H0 = sys.H0
for i in 1:2^L
    for j in 1:2^L
        val = MatrixChecker(H0, sys.sites, i-1, j-1) # -1 because of 0-based indexing in binary representation
        if abs(val) > 1e-6
            print(@sprintf "<%.0f|H0|%.0f> = %8.3f   " i-1 j-1 val)
        end
    end
    println()
end



# 3. Purification Check


