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
include("../src/hf/self_consistent.jl")

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)


println("\n--- Test 1: System Construction ---")
L = 6
τ = (sqrt(10)-2.0)/2.0
t(x) = -1.0 - 0.1 * cos(2π * τ * (x-0.5))
U = 10.0
W = nothing
S(x) = 0.5 * cos(π * x)

tci_tol = 1e-6
itensors_tol = 1e-10
itensors_maxdim = 100
density = 0.5
purification_steps = 100
scf_mixing = 0.9
scf_tol = 0.1 #%
scf_max_iterations = 100

params = ModelParameters(L, t, U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys = System(params)

println("System constructed successfully:")

show(sys)

H_min = -20.0
H_max = 20.0

run_scf!(sys,H_min,H_max, verbose=:all)


println("\n--- Final Hartree Potential Check ---")
print("\nFinal Hartree MPO values...\n")
for i in 1:2^L
    val = MatrixChecker(sys.ρ, sys.sites, i , i , sys.bra_states, sys.ket_states)
    if abs(val) > 1e-6
        println(@sprintf "<%.0f|ρ|%.0f> = %8.3f   " i-1 i-1 val)
    end
end
println()   
