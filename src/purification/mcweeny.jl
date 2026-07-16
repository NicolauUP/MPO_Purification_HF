function construct_rho_0_mcweeny(sys::System, params::AbstractModelParameters, mu, H_min::Float64, H_max::Float64; to_gpu=identity)
    Id = to_gpu(Identity_MPO(sys.sites))
    c = 0.5 / max(mu - H_min, H_max - mu)
    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return +((0.5 + c * mu) * Id, -c * H;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

function perform_purification_mcweeny_mu(rho::MPO, params::AbstractModelParameters; verbose::Int=1, spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing, spectral_bounds_validation::Symbol=:not_provided, chemical_potential::Union{Nothing,Real}=nothing)
    isnothing(chemical_potential) && throw(ArgumentError("method=:mcweeny_mu requires chemical_potential"))
    isnothing(spectral_bounds) && throw(ArgumentError("method=:mcweeny_mu requires spectral_bounds"))
    bounds = validate_spectral_bounds(spectral_bounds...)
    bounds[1] < chemical_potential < bounds[2] || throw(ArgumentError("chemical_potential must lie inside spectral_bounds"))
    for iteration in 1:params.purification_steps
        P2 = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        idem = idempotency_residual(rho, P2)
        if idem < 1e-3
            return purification_result(rho, params; method=:mcweeny_mu, converged=true, termination_reason=:idempotency_threshold, iterations=iteration, spectral_bounds=bounds, spectral_bounds_validation=spectral_bounds_validation, target_particles=nothing)
        end
        P3 = apply(rho, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        rho = +(3.0 * P2, -2.0 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end
    @warn "McWeeny purification did not converge after $(params.purification_steps) steps."
    return purification_result(rho, params; method=:mcweeny_mu, converged=false, termination_reason=:max_iterations, iterations=params.purification_steps, spectral_bounds=bounds, spectral_bounds_validation=spectral_bounds_validation, target_particles=nothing)
end
