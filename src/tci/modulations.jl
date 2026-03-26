"""
    build_pi_modulation(sys, W0; ϵ)

Build a diagonal MPO for a π-modulated on-site potential:
    W(i) = W0 * cos(π * i)
defined on integer sites i ∈ {0, ..., 2^L - 1}.
"""
function build_pi_modulation(sys::System, W0::Float64; ϵ::Float64=1e-10)
    f = x -> W0 * cos( π * x)
    _, mpo, _ = Quantics_TCI(f, Float64, sys.sites, ϵ)
    return mpo
end