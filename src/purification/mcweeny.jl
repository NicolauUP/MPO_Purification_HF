function construct_rho_0_mcweeny(sys::System, params::AbstractModelParameters, mu, H_min::Float64, H_max::Float64; to_gpu=identity)
    Id = to_gpu(Identity_MPO(sys.sites))
    c = 0.5 / max(mu - H_min, H_max - mu)
    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return +((0.5 + c * mu) * Id, -c * H;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

"""Apply one McWeeny update to an already computed square.

`:standard` reproduces `3P²-2P³`. `:horner` evaluates the identical
polynomial as `P²(3I-2P)`, avoiding a separately compressed `P³`. The Horner
form is opt-in because finite-rank MPO truncation makes the two evaluation
orders numerically distinct.
"""
function _mcweeny_update_from_square(
    rho::MPO,
    rho_squared::MPO,
    params::AbstractModelParameters;
    form::Symbol=:standard,
    identity_mpo::Union{Nothing,MPO}=nothing,
)
    if form == :standard
        rho_cubed = apply(rho, rho_squared;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        rho_next = +(3.0 * rho_squared, -2.0 * rho_cubed;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        return (
            rho=rho_next,
            auxiliary=rho_cubed,
            auxiliary_kind=:rho_cubed,
        )
    elseif form == :horner
        isnothing(identity_mpo) && throw(ArgumentError(
            "form=:horner requires identity_mpo",
        ))
        factor = +(3.0 * identity_mpo, -2.0 * rho;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        rho_next = apply(rho_squared, factor;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        return (
            rho=rho_next,
            auxiliary=factor,
            auxiliary_kind=:horner_factor,
        )
    end
    throw(ArgumentError("unknown McWeeny polynomial form $form; use :standard or :horner"))
end

function perform_purification_mcweeny_mu(rho::MPO, params::AbstractModelParameters; verbose::Int=1, spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing, spectral_bounds_validation::Symbol=:not_provided, chemical_potential::Union{Nothing,Real}=nothing)
    isnothing(chemical_potential) && throw(ArgumentError("method=:mcweeny_mu requires chemical_potential"))
    isnothing(spectral_bounds) && throw(ArgumentError("method=:mcweeny_mu requires spectral_bounds"))
    bounds = validate_spectral_bounds(spectral_bounds...)
    bounds[1] < chemical_potential < bounds[2] || throw(ArgumentError("chemical_potential must lie inside spectral_bounds"))
    # `:mcweeny_mu` has no prescribed particle number to check independently.
    # A 1e-3 idempotency residual can still leave O(1e-5) occupation errors on
    # small finite systems, so require the next cubic-convergence regime.
    idempotency_tolerance = 1e-6
    for iteration in 1:params.purification_steps
        P2 = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        idem = idempotency_residual(rho, P2)
        if idem < idempotency_tolerance
            return purification_result(rho, params; method=:mcweeny_mu, converged=true, termination_reason=:idempotency_threshold, iterations=iteration, spectral_bounds=bounds, spectral_bounds_validation=spectral_bounds_validation, target_particles=nothing)
        end
        rho = _mcweeny_update_from_square(rho, P2, params).rho
    end
    @warn "McWeeny purification did not converge after $(params.purification_steps) steps."
    return purification_result(rho, params; method=:mcweeny_mu, converged=false, termination_reason=:max_iterations, iterations=params.purification_steps, spectral_bounds=bounds, spectral_bounds_validation=spectral_bounds_validation, target_particles=nothing)
end
