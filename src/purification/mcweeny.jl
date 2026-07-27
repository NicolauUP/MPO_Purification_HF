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

function _validate_mcweeny_trace_guard(
    trace_target::Union{Nothing,Real},
    trace_tolerance::Union{Nothing,Real},
)
    isnothing(trace_target) == isnothing(trace_tolerance) || throw(ArgumentError(
        "trace_target and trace_tolerance must be supplied together",
    ))
    if !isnothing(trace_target)
        isfinite(trace_target) && trace_target >= 0 || throw(ArgumentError(
            "trace_target must be finite and nonnegative",
        ))
        isfinite(trace_tolerance) && trace_tolerance >= 0 || throw(ArgumentError(
            "trace_tolerance must be finite and nonnegative",
        ))
    end
    return nothing
end

"""Purify a fixed-chemical-potential density matrix with the McWeeny map.

Every evaluated iterate is ranked by idempotency. If the production
`1e-6` threshold is not reached, the best evaluated iterate is returned
instead of the final, potentially deteriorated, finite-rank MPO.

`trace_target` and `trace_tolerance` are an optional eligibility guard and
must be supplied together. They are never inferred from `params.density`,
since doing so would change this fixed-chemical-potential method into a
canonical one. `form=:horner` additionally requires an identity MPO on the
same device as `rho`.
"""
function perform_purification_mcweeny_mu(
    rho::MPO,
    params::AbstractModelParameters;
    verbose::Int=1,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
    chemical_potential::Union{Nothing,Real}=nothing,
    form::Symbol=:standard,
    identity_mpo::Union{Nothing,MPO}=nothing,
    trace_target::Union{Nothing,Real}=nothing,
    trace_tolerance::Union{Nothing,Real}=nothing,
)
    isnothing(chemical_potential) && throw(ArgumentError("method=:mcweeny_mu requires chemical_potential"))
    isnothing(spectral_bounds) && throw(ArgumentError("method=:mcweeny_mu requires spectral_bounds"))
    bounds = validate_spectral_bounds(spectral_bounds...)
    bounds[1] < chemical_potential < bounds[2] || throw(ArgumentError("chemical_potential must lie inside spectral_bounds"))
    form in (:standard, :horner) || throw(ArgumentError(
        "form must be :standard or :horner",
    ))
    form == :horner && isnothing(identity_mpo) && throw(ArgumentError(
        "form=:horner requires identity_mpo",
    ))
    _validate_mcweeny_trace_guard(trace_target, trace_tolerance)

    # `:mcweeny_mu` has no prescribed particle number to check independently.
    # A 1e-3 idempotency residual can still leave O(1e-5) occupation errors on
    # small finite systems, so require the next cubic-convergence regime.
    idempotency_tolerance = 1e-6
    best_rho = nothing
    best_idempotency = Inf
    best_iteration = 0
    last_evaluated_rho = rho
    max_bond_dimension = maxlinkdim(rho)
    bond_dimension_sum = 0
    bond_dimension_samples = 0
    bond_dimensions = NTuple{3,Int}[]
    cubes = 0

    for iteration in 1:params.purification_steps
        P2 = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        idem = idempotency_residual(rho, P2)
        trace_value = real(tr(rho))
        trace_eligible = isnothing(trace_target) ||
            abs(trace_value - trace_target) <= trace_tolerance
        if trace_eligible && isfinite(idem) && idem < best_idempotency
            best_rho = rho
            best_idempotency = idem
            best_iteration = iteration
        end
        last_evaluated_rho = rho

        auxiliary = nothing
        if !(idem < idempotency_tolerance && trace_eligible) &&
           iteration < params.purification_steps
            update = _mcweeny_update_from_square(
                rho, P2, params;
                form=form,
                identity_mpo=identity_mpo,
            )
            auxiliary = update.auxiliary
            rho = update.rho
            cubes += form == :standard
        end

        rho_chi = maxlinkdim(last_evaluated_rho)
        square_chi = maxlinkdim(P2)
        auxiliary_chi = isnothing(auxiliary) ? 0 : maxlinkdim(auxiliary)
        max_bond_dimension = max(
            max_bond_dimension, rho_chi, square_chi, auxiliary_chi,
        )
        bond_dimension_sum += rho_chi + square_chi
        bond_dimension_samples += 2
        if !isnothing(auxiliary)
            bond_dimension_sum += auxiliary_chi
            bond_dimension_samples += 1
        end
        push!(bond_dimensions, (rho_chi, square_chi, auxiliary_chi))
        work = PurificationWorkStats(
            iteration,
            cubes,
            max_bond_dimension,
            bond_dimension_sum / bond_dimension_samples,
            copy(bond_dimensions),
        )

        if idem < idempotency_tolerance && trace_eligible
            return purification_result(
                last_evaluated_rho, params;
                method=:mcweeny_mu,
                converged=true,
                termination_reason=:idempotency_threshold,
                iterations=iteration,
                selected_iteration=iteration,
                spectral_bounds=bounds,
                spectral_bounds_validation=spectral_bounds_validation,
                target_particles=trace_target,
                work=work,
            )
        end
    end

    selected_rho = isnothing(best_rho) ? last_evaluated_rho : best_rho
    selected_iteration = isnothing(best_rho) ? params.purification_steps : best_iteration
    termination_reason = if isnothing(best_rho)
        :trace_guard_rejected_all
    elseif best_iteration < params.purification_steps
        :best_iterate_after_max_iterations
    else
        :max_iterations
    end
    @warn "McWeeny purification did not converge after $(params.purification_steps) steps; returning selected iteration $selected_iteration."
    return purification_result(
        selected_rho, params;
        method=:mcweeny_mu,
        converged=false,
        termination_reason=termination_reason,
        iterations=params.purification_steps,
        selected_iteration=selected_iteration,
        spectral_bounds=bounds,
        spectral_bounds_validation=spectral_bounds_validation,
        target_particles=trace_target,
        work=PurificationWorkStats(
            params.purification_steps,
            cubes,
            max_bond_dimension,
            bond_dimension_sum / max(1, bond_dimension_samples),
            copy(bond_dimensions),
        ),
    )
end
