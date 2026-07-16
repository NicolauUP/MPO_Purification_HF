function construct_rho_0_mcweeny(sys::System, params::AbstractModelParameters, mu, H_min::Float64, H_max::Float64; to_gpu=identity)
    Id = to_gpu(Identity_MPO(sys.sites))
    c = 0.5 / max(mu - H_min, H_max - mu)
    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return +((0.5 + c * mu) * Id, -c * H;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

"""Legacy fixed-particle-number outer-μ search around standard McWeeny."""
function perform_purification_grandcanonical(sys::System, params::AbstractModelParameters, H_min::Float64, H_max::Float64; verbose::Int=1, to_gpu=identity)
    Ne = round(Int, 2^params.L * params.density)
    mu_low, mu_high = -1.0, 1.0
    rho_new = nothing
    while (mu_high - mu_low) > 1e-4
        mu = (mu_high + mu_low) / 2.0
        rho_0 = construct_rho_0_mcweeny(sys, params, mu, H_min, H_max; to_gpu=to_gpu)
        for _ in 1:params.purification_steps
            P2 = apply(rho_0, rho_0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            P3 = apply(rho_0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            rho_new = +(3.0 * P2, -2.0 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            idempotency_residual(rho_0, P2) < 1e-3 && break
            rho_0 = rho_new
        end
        trace_value = real(tr(rho_new))
        abs(trace_value - Ne) < 0.5 && return rho_new
        trace_value > Ne ? (mu_high = mu) : (mu_low = mu)
    end
    @warn "Legacy grandcanonical μ-search reached its resolution limit."
    return rho_new
end
