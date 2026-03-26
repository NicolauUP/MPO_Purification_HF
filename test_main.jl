# test_main.jl
__precompile__(false)
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
import TensorCrossInterpolation as TCI
using LinearAlgebra

include("src/core/operators.jl")
include("src/core/system.jl")
# include("src/hamiltonians/mpo_construction.jl")
# include("src/purification/mcweeny.jl")
# include("src/utils/quantics.jl")
# include("src/tci/modulations.jl")

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)

# 1. System construction
# ─────────────────────────────────────────
println("\n--- Test 1: System Construction ---")
L = 4
t = 1.0
U = 0.0
W = nothing
params = ModelParameters(L, t, U, W)
sys = System(params)
println("System constructed successfully:")
println()
show(sys)

