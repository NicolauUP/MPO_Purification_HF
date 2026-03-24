
"""
    construct_rho_0(sys, H, ϵ, maxχ, H_max, H_min, Ne)

Build the initial density matrix guess by linearly mapping the
eigenvalues of H into [0,1] with the correct electron count Ne.
"""
function construct_rho_0(sys::System, H::MPO, ϵ::Float64, maxχ::Int,
                          H_max::Float64, H_min::Float64, Ne::Int)
    N = length(sys.sites)
    Id = Identity_MPO(sys.sites)   # built internally!
    μ = real(tr(H) / N)
    λ = minimum((Ne / (H_max - μ), (N - Ne) / (μ - H_min)))
    coeff_I = (Ne + λ * μ) / N
    coeff_H = -(λ / N)
    return +(coeff_I * Id, coeff_H * H; cutoff=ϵ, maxdim=maxχ)
end

