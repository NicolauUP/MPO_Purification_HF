# test_pi_modulation.jl
using ITensors, ITensorMPS
using QuanticsTCI, TensorCrossInterpolation
import TensorCrossInterpolation as TCI
import QuanticsTCI: Quantics

include("src/core/operators.jl")
include("src/core/system.jl")
include("src/utils/quantics.jl")
include("src/tci/modulations.jl")

# System parameters
L  = 8
n  = 1
t  = 1.0
U  = 0.0
W  = nothing

sys = System(1, n, t, U, W, L)

# Build π-modulation MPO
W0  = 1.0
mpo = build_pi_modulation(sys, W0; ϵ=1e-10)

# Verify diagonal elements against exact values
println("Verifying π-modulation MPO diagonal elements:")
println("─────────────────────────────────────────────")
for i in 0:2^L-1
    val_mpo   = real(MatrixChecker(mpo, sys.sites, i, i))
    val_exact = W0 * cos(π * i)
    err       = abs(val_mpo - val_exact)
    println("i=$i | MPO: $val_mpo | Exact: $val_exact | Error: $err")
end