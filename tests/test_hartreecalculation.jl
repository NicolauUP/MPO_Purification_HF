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
include("../src/purification/mcweeny.jl")
include("../src/tci/modulations.jl")
include("../src/tci/density_matrix.jl")


println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)


println("\n--- Test 1: System Construction ---")
L = 6
t = -1.0
U = 1.0
W = nothing
S(x) = 0.5 * cos(pi * x)

tci_tol = 1e-6
itensors_tol = 1e-10
itensors_maxdim = 100
density = 0.5
purification_steps = 40

params = ModelParameters(L, t, U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps)
sys = System(params)

println("System constructed successfully:")
println()
show(sys)

ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
println("\n--- Test 3: Purification Check ---")
ρ_purified = perform_purification(ρ0, params; verbose=1)

sys.ρ = ρ_purified
println("\n--- Test 2: MPO Check ---")
@time vh_mpo = extract_hartree_mpo_1d(sys)

print("\nHartree MPO extracted successfully. Checking values...\n")
for i in 1:2^L
    val = MatrixChecker(vh_mpo, sys.sites, i , i , sys.bra_states, sys.ket_states)
    if abs(val) > 1e-6
        println(@sprintf "<%.0f|VH|%.0f> = %8.3f   " i-1 i-1 val)
    end
end
println()
