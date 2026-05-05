


using MPO_MeanField
using CUDA
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
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
println("HartreeFockMPO — Timings in 2D")
println("="^50)


τ = (sqrt(10)-2.0)/2.0





L = parse(Int, ARGS[1])
τ = (sqrt(10)-2.0)/2.0
V2 = parse(Float64, ARGS[2])
W_Amp = parse(Float64, ARGS[3])
maxdim = parse(Int, ARGS[4])
tol = parse(Float64, ARGS[5])
function tx_qp(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return -1.0 - V2 * (cos(2*pi * τ * (x-0.5)))
end
function ty_qp(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return -1.0 - V2 * (cos(2*pi * τ * (y-0.5)))
end


U = 0.3   
W = nothing
function S(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return W_Amp * (cos(pi*y) * cos(pi*x))
end


tci_tol = 1e-6
itensors_tol = tol
itensors_maxdim = maxdim
density = 0.5
purification_steps = 200
scf_mixing = 0.9
scf_tol = 0.1 #%
scf_max_iterations = 100

params = ParametersSquare(L, (tx_qp,ty_qp), U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys = System(params)
show(sys)


# ─────────────────────────────────────────
# JIT Warmup — compile all functions on a small system
# ─────────────────────────────────────────
println("\n--- JIT Warmup (L=4) ---")
params_warmup = ParametersSquare(4, (tx_qp,ty_qp), U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys_warmup = System(params_warmup)

ρ0_warmup = construct_rho_0(sys_warmup, params_warmup, -3.0, 3.0)
ρ0_warmup = to_gpu(ρ0_warmup)
_ = perform_purification(ρ0_warmup, params_warmup; verbose=0)
ρ0_warmup = nothing
GC.gc()
CUDA.functional() && CUDA.reclaim()
println("Warmup complete. JIT compilation done.\n")

# ─────────────────────────────────────────
# Timed real run — L=20
# ─────────────────────────────────────────
L = 30
params = ParametersSquare(L, (tx_qp,ty_qp), U, W,S, tci_tol, itensors_tol, itensors_maxdim, density, purification_steps, scf_mixing, scf_tol, scf_max_iterations)
sys = System(params)
println("\n--- Test 1a: Initial ρ₀ Trace ---")
t_rho0 = @elapsed begin
    ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
end
T1 = real(tr(ρ0))
Ne = round(Int, 2^L * density)
@assert isapprox(T1, Ne) "Trace wrong: got $T1, expected $Ne"
println("Tr(ρ₀) = $T1 ✓  [$(round(t_rho0, digits=2))s]")

println("\n--- Test 2a: Purified ρ Idempotency ---")
ρ0 = to_gpu(ρ0)
t_purification = @elapsed begin
    ρ_purified = perform_purification(ρ0, params; verbose=1)
end
ρ_purified = to_cpu(ρ_purified)
println("Purification time: $(round(t_purification, digits=2))s")

# ... rest of tests with timing
println("\n--- Timing Summary ---")
println(@sprintf "  construct_rho_0    : %8.2f s" t_rho0)
println(@sprintf "  perform_purification: %8.2f s" t_purification)

