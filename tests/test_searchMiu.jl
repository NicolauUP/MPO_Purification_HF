using MPO_MeanField
using CUDA
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf

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

U = 0.3
W = nothing
S(x) = 0.25 * cos(π * x)
 
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
L = 20

# Corrected to use Parameters1D and keyword arguments
params = Parameters1D(
    L=L, t=t, U=U, W=W, S=S, 
    tci_tol=tci_tol, itensors_tol=itensors_tol, 
    itensors_maxdim=itensors_maxdim, density=density, 
    purification_steps=purification_steps, scf_mixing=scf_mixing, 
    scf_tol=scf_tol, scf_max_iterations=scf_max_iterations
)

sys = System(params)




t_purification = @elapsed begin
    ρ_purified = perform_purification_grandcanonical(sys, params, -5.0, 5.0; verbose=1, to_gpu=to_gpu)
end
ρ_purified = to_cpu(ρ_purified)
println("Purification time: $(round(t_purification, digits=2))s")

# ─────────────────────────────────────────
# Timing Summary
# ─────────────────────────────────────────
println("\n--- Timing Summary ---")
println(@sprintf "  perform_purification : %8.2f s" t_purification)