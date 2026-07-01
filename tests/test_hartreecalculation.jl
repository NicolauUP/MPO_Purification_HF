using MPO_MeanField
using CUDA
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf
using HDF5

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
println("HartreeFockMPO — Search of the chemical potential in McWeeny purification")
println("="^50)

τ = (sqrt(10)-2.0)/2.0
t(x) = -1.0 - 2.0 * cos(2π * τ * (x - 0.5))
#t = -1.0
U = 0.3
W = nothing
S(x) = 0.1 * cos(π * x)
 
tci_tol = 1e-6
itensors_tol = 1e-10
itensors_maxdim = 600
density = 0.5
purification_steps = 200
scf_mixing = 0.9
scf_tol = 0.1 # %
scf_max_iterations = 100

# ─────────────────────────────────────────
# JIT Warmup — compile all functions on a small system
# ─────────────────────────────────────────
println("\n--- JIT Warmup (L=4) ---")

# Corrected to use Parameters1D and keyword arguments
params_warmup = Parameters1D(
    L=4, t=t, U=U, W=W, S=S, 
    tci_tol=tci_tol, itensors_tol=itensors_tol, 
    itensors_maxdim=itensors_maxdim, density=density, 
    purification_steps=10, scf_mixing=scf_mixing, 
    scf_tol=scf_tol, scf_max_iterations=scf_max_iterations
)

sys_warmup = System(params_warmup)
rho = perform_purification_grandcanonical(sys_warmup, params_warmup, -5.0,5.0; verbose=1, to_gpu=to_gpu)

rho = nothing
GC.gc()
CUDA.functional() && CUDA.reclaim()
println("Warmup complete. JIT compilation done.\n")

# ─────────────────────────────────────────
# Timed real run — L=30
# ─────────────────────────────────────────
L = 10

# Corrected to use Parameters1D and keyword arguments
params = Parameters1D(
    L=L, t=t, U=U, W=W, S=S, 
    tci_tol=tci_tol, itensors_tol=itensors_tol, 
    itensors_maxdim=itensors_maxdim, density=density, 
    purification_steps=purification_steps, scf_mixing=scf_mixing, 
    scf_tol=scf_tol, scf_max_iterations=scf_max_iterations
)

sys = System(params)


# 1) Compute density matrix - CHECK
rho = perform_purification_grandcanonical(sys, params, -5.0,5.0; verbose=1, to_gpu=to_gpu) 

# 2.1) Compute Hartree potential for the dense matrix of rho
ni_matrix = zeros(ComplexF64, 2^L)
for i in 1:2^L
    ni_matrix[i] = MatrixChecker(rho, sys.sites, i,i, sys.bra_states, sys.ket_states)
end

H_potential_dense = zeros(ComplexF64, 2^L)
H_potential_dense[1] = U * ni_matrix[2]
H_potential_dense[end] = U * ni_matrix[end-1]
for i in 2:2^L-1
    H_potential_dense[i] = U * (ni_matrix[i-1] + ni_matrix[i+1])
end

sys.ρ = rho

output_file = "Tests_HartreeFock_MPO_L$(L).h5"
h5open(output_file, "w") do file
    write(file, "H_potential_dense", H_potential_dense)
    write(file, "ni_matrix", ni_matrix)
end

hartree_mpo = extract_hartree_mpo_1d(sys)

H_potential_mpo = zeros(ComplexF64, 2^L)
for i in 1:2^L
    H_potential_mpo[i] = MatrixChecker(hartree_mpo, sys.sites, i,i, sys.bra_states, sys.ket_states)
end 
h5open(output_file, "r+") do file
    write(file, "H_potential_mpo", H_potential_mpo)
end

differences = abs.(H_potential_dense - H_potential_mpo)
h5open(output_file, "r+") do file
    write(file, "differences_HMPO_twoways", differences) 
end



# TODO:
#= 
 1) Compute density matrix - CHECK
 2.1) Compute Hartree potencial for the dense matrix of rho
 2.2) Compute Hartree potential from the contractions and Quantics
 2.3) Compare the two Hartree potentials

 3.1) Compute Fock potencial for the dense matrix of rho
 3.2) Compute Fock potential from the contractions and Quantics
 3.3) Compare the two Fock potentials

 4) Verify Mean Field hamiltonian MPO it's the same!

 5) Run self-consistent loop and check convergence!

=#