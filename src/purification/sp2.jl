function _sp2_trace_tolerance(params::AbstractModelParameters, Ne::Int)
    return max(10sqrt(eps(Float64)), 10params.itensors_tol) * max(1, Ne)
end

function _sp2_hermiticity_tolerance(params::AbstractModelParameters)
    return max(1e-10, 10params.itensors_tol)
end

function _sp2_branch(trace_value::Real, trace_squared::Real, Ne::Int, trace_tolerance::Real)
    trace_hole = 2trace_value - trace_squared
    if trace_value > Ne + trace_tolerance
        return :square
    elseif trace_value < Ne - trace_tolerance
        return :hole
    elseif abs(trace_squared - Ne) <= abs(trace_hole - Ne) + trace_tolerance
        return :square
    end
    return :hole
end

"""
    perform_purification_sp2(rho0, params; ...)

Canonical trace-correcting second-order purification. Each iteration forms one
compressed square and selects either `X²` or `2X-X²` by its trace relative to
the requested occupation.
"""
function perform_purification_sp2(
    rho0::MPO,
    params::AbstractModelParameters;
    verbose::Int=1,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    idempotency_tolerance::Real=1e-3,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
)
    !isnothing(spectral_bounds) || throw(ArgumentError(
        "SP2 requires explicit enclosing spectral_bounds",
    ))
    spectral_bounds = validate_spectral_bounds(spectral_bounds...)
    idempotency_tolerance = Float64(idempotency_tolerance)
    isfinite(idempotency_tolerance) && idempotency_tolerance > 0 ||
        throw(ArgumentError("SP2 idempotency_tolerance must be finite and positive"))

    N = 2^params.L
    Ne = round(Int, N * params.density)
    0 < Ne < N || throw(ArgumentError("SP2 supports only 0 < Ne < N, got Ne=$Ne"))
    trace_tolerance = _sp2_trace_tolerance(params, Ne)
    hermiticity_tolerance = _sp2_hermiticity_tolerance(params)

    verbose > 0 && println(io, "SP2 purifying N=$N, density=$(params.density), Ne=$Ne")
    rho = rho0
    previous_idempotency = Inf
    stagnant_steps = 0
    max_bond_dimension = maxlinkdim(rho)
    bond_dimension_sum = 0
    bond_dimension_samples = 0
    bond_dimensions = NTuple{3,Int}[]

    for iteration in 1:params.purification_steps
        rho_squared = apply(rho, rho;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        max_bond_dimension = max(max_bond_dimension, maxlinkdim(rho), maxlinkdim(rho_squared))
        bond_dimension_sum += maxlinkdim(rho) + maxlinkdim(rho_squared)
        bond_dimension_samples += 2
        push!(bond_dimensions, (maxlinkdim(rho), maxlinkdim(rho_squared), 0))
        work = PurificationWorkStats(
            iteration, 0, max_bond_dimension, bond_dimension_sum / bond_dimension_samples,
            copy(bond_dimensions),
        )
        trace_value = real(tr(rho))
        trace_squared = real(tr(rho_squared))
        idem = idempotency_residual(rho, rho_squared)
        herm = _relative_mpo_residual(rho, ITensors.dag(rho), params)

        if verbose > 0
            details = @sprintf(
                "Tr=%.12g | err=%.3e | idem=%.3e | herm=%.3e | χ=(%d,%d)",
                trace_value,
                abs(trace_value - Ne),
                idem,
                herm,
                maxlinkdim(rho),
                maxlinkdim(rho_squared),
            )
            print_iteration_progress(
                io, "SP2", iteration, params.purification_steps, details;
                overwrite=overwrite_progress,
            )
        end

        if abs(trace_value - Ne) <= trace_tolerance &&
           idem < idempotency_tolerance && herm <= hermiticity_tolerance
            verbose > 0 && finish_iteration_progress(io, overwrite_progress)
            return purification_result(
                rho, params;
                method=:sp2,
                converged=true,
                termination_reason=:idempotency_threshold,
                iterations=iteration,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
                work=work,
            )
        end

        if idem >= previous_idempotency * (1 - 1e-12)
            stagnant_steps += 1
        else
            stagnant_steps = 0
        end
        if stagnant_steps >= 8
            verbose > 0 && finish_iteration_progress(io, overwrite_progress)
            return purification_result(
                rho, params;
                method=:sp2,
                converged=false,
                termination_reason=:stagnation,
                iterations=iteration,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
                work=work,
            )
        end
        previous_idempotency = idem

        branch = _sp2_branch(trace_value, trace_squared, Ne, trace_tolerance)
        if branch == :square
            rho = rho_squared
        else
            rho = +(2.0 * rho, -rho_squared;
                cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
            )
        end
    end

    verbose > 0 && finish_iteration_progress(io, overwrite_progress)
    @warn "SP2 purification did not converge after $(params.purification_steps) steps."
    return purification_result(
        rho, params;
        method=:sp2,
        converged=false,
        termination_reason=:max_iterations,
        iterations=params.purification_steps,
        spectral_bounds=spectral_bounds,
        spectral_bounds_validation=spectral_bounds_validation,
        work=PurificationWorkStats(
            params.purification_steps, 0, max_bond_dimension,
            bond_dimension_sum / max(1, bond_dimension_samples),
            copy(bond_dimensions),
        ),
    )
end
