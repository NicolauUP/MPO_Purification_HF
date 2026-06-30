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

params = Parameters1D(L, t, U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps)
sys = System(params)


# TODO:
#= 
 1) Compute density matrix 
 2.1) Compute Hartree potencial for the dense matrix of rho
 2.2) Compute Hartree potential from the contractions and Quantics
 2.3) Compare the two Hartree potentials

 3.1) Compute Fock potencial for the dense matrix of rho
 3.2) Compute Fock potential from the contractions and Quantics
 3.3) Compare the two Fock potentials

 4) Verify Mean Field hamiltonian MPO it's the same!

 5) Run self-consistent loop and check convergence!

=#