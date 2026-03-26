# test_main.jl
__precompile__(false)
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using LinearAlgebra

include("src/core/operators.jl")
include("src/core/system.jl")
include("src/hamiltonians/mpo_construction.jl")
include("src/purification/mcweeny.jl")
include("src/utils/quantics.jl")
include("src/tci/modulations.jl")

println("="^50)
println("HartreeFockMPO — Test Suite")
println("="^50)

# ─────────────────────────────────────────
# Parameters
# ─────────────────────────────────────────
L   = 4
n   = 1
d   = 1
t   = 1.0
U   = 0.0
W   = nothing
Ne  = L ÷ 2      # half filling
H_max =  2.0 * t + 5.0
H_min = -2.0 * t - 5.0

# ─────────────────────────────────────────
# 1. System construction
# ─────────────────────────────────────────
println("\n--- Test 1: System Construction ---")
sys = System(d, n, t, U, W, L)
println(sys)

# ─────────────────────────────────────────
# 2. H0 construction
# ─────────────────────────────────────────
println("\n--- Test 2: H0 Construction ---")
H0 = build_H0_chain(sys; cutoff=1e-10, maxdim=100)
println("H0 built successfully")
println("Bond dimension: ", maxlinkdim(H0))

# ─────────────────────────────────────────
# 3. H0 matrix elements
# ─────────────────────────────────────────
println("\n--- Test 3: H0 Matrix Elements ---")
println("Full $(2^L) x $(2^L) matrix:")
for i in 0:2^L-1
    for j in 0:2^L-1
        val = real(MatrixChecker(H0, sys.sites, i, j))
        if abs(val) > 1e-10
            print("  ⟨$i|H0|$j⟩ = $(round(val,digits=2))")
        end
    end
    println()
end

# ─────────────────────────────────────────
# 4. π-modulation MPO
# ─────────────────────────────────────────
println("\n--- Test 4: π-modulation MPO ---")
W0      = 10.0
W_mpo   = build_pi_modulation(sys, W0; ϵ=1e-10)
println("π-modulation MPO built successfully")
println("Bond dimension: ", maxlinkdim(W_mpo))

println("\nDiagonal elements vs exact cos(π*i):")
max_err = 0.0
for i in 0:2^L-1
    val_mpo   = real(MatrixChecker(W_mpo, sys.sites, i, i))
    val_exact = W0 * cos(π  * i)
    err       = abs(val_mpo - val_exact)
    global max_err = max(max_err, err)
    println("  i=$i | MPO: $(round(val_mpo, digits=8)) | Exact: $(round(val_exact, digits=8)) | Error: $(round(err, digits=2))")
end
println("Max error: $max_err")

# ─────────────────────────────────────────
# 5. Purification
# ─────────────────────────────────────────
println("\n--- Test 5: McWeeny Purification ---")
ρ0 = construct_rho_0(sys, H0, 1e-10, 100, H_max, H_min, Ne)
println("ρ0 built, initial trace: ", real(tr(ρ0)))

ρ = perform_purification(ρ0; maxχ=100, ϵ=1e-10, max_steps=40, verbose=1)

println("\nFinal trace (should be Ne=$Ne): ", real(tr(ρ)))
println("Bond dimension of ρ: ", maxlinkdim(ρ))

# ─────────────────────────────────────────
# 6. Density matrix diagonal
# ─────────────────────────────────────────
println("\n--- Test 6: Density Matrix Diagonal ---")
println("Occupation per basis state ⟨i|ρ|i⟩:")
for i in 0:2^L-1
    val = real(MatrixChecker(ρ, sys.sites, i, i))
    if abs(val) > 1e-10
        println("  ⟨$i|ρ|$i⟩ = $(round(val, digits=8))")
    end
end

println("\n", "="^50)
println("All tests complete!")
println("="^50)