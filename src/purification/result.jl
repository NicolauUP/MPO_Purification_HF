"""
    PurificationResult

Diagnostics returned by a density-matrix purification backend. `spectral_bounds`
is `nothing` when the backend was not given explicit bounds.
"""
struct PurificationResult{R}
    rho::R
    method::Symbol
    converged::Bool
    termination_reason::Symbol
    iterations::Int
    trace::Float64
    trace_error::Float64
    idempotency_residual::Float64
    hermiticity_residual::Float64
    final_bond_dimension::Int
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}
    spectral_bounds_validation::Symbol
    truncation_cutoff::Float64
    maxdim::Int
end

function validate_spectral_bounds(H_min::Real, H_max::Real; safety_margin::Real=0.0)
    isfinite(H_min) || throw(ArgumentError("H_min must be finite, got $H_min"))
    isfinite(H_max) || throw(ArgumentError("H_max must be finite, got $H_max"))
    H_min < H_max || throw(ArgumentError("spectral bounds must satisfy H_min < H_max, got [$H_min, $H_max]"))
    isfinite(safety_margin) && safety_margin >= 0 ||
        throw(ArgumentError("safety_margin must be finite and nonnegative, got $safety_margin"))
    return (Float64(H_min), Float64(H_max))
end

"""
    verify_spectral_bounds_exact(sys, H, H_min, H_max; safety_margin=0.0)

Validate user-supplied bounds against the dense spectrum of `H`. This is a
CPU-only small-system test oracle and intentionally rejects systems with more
than 16 basis states; it is not a production spectral estimator.
"""
function verify_spectral_bounds_exact(
    sys::System,
    H::MPO,
    H_min::Real,
    H_max::Real;
    safety_margin::Real=0.0,
)
    bounds = validate_spectral_bounds(H_min, H_max; safety_margin=safety_margin)
    N = 2^length(sys.sites)
    N <= 16 || throw(ArgumentError(
        "exact spectral-bound validation is limited to N <= 16; got N=$N",
    ))

    matrix = [
        MatrixChecker(H, sys.sites, i, j, sys.bra_states, sys.ket_states)
        for i in 1:N, j in 1:N
    ]
    hermiticity_error = opnorm(matrix - matrix')
    tolerance = sqrt(eps(Float64)) * max(1.0, opnorm(matrix))
    hermiticity_error <= tolerance || throw(ArgumentError(
        "cannot validate spectral bounds for a non-Hermitian effective Hamiltonian " *
        "(residual=$hermiticity_error)",
    ))

    eigenvalues = eigvals(Hermitian((matrix + matrix') / 2))
    λ_min, λ_max = extrema(real.(eigenvalues))
    lower_required = λ_min - safety_margin
    upper_required = λ_max + safety_margin
    bounds[1] <= lower_required + tolerance || throw(ArgumentError(
        "H_min=$(bounds[1]) does not enclose λ_min=$λ_min with safety_margin=$safety_margin",
    ))
    bounds[2] >= upper_required - tolerance || throw(ArgumentError(
        "H_max=$(bounds[2]) does not enclose λ_max=$λ_max with safety_margin=$safety_margin",
    ))
    return (Float64(λ_min), Float64(λ_max))
end

function _relative_mpo_residual(A::MPO, B::MPO, params::AbstractModelParameters)
    difference = +(A, -B; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    numerator = sqrt(max(0.0, real(inner(difference, difference))))
    denominator = max(sqrt(max(0.0, real(inner(B, B)))), sqrt(eps(Float64)))
    return numerator / denominator
end

function purification_result(
    rho::MPO,
    params::AbstractModelParameters;
    method::Symbol,
    converged::Bool,
    termination_reason::Symbol,
    iterations::Int,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
)
    Ne = round(Int, 2^params.L * params.density)
    rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    trace_value = Float64(real(tr(rho)))
    idempotency_residual = _relative_mpo_residual(rho_squared, rho, params)
    hermiticity_residual = _relative_mpo_residual(rho, ITensors.dag(rho), params)

    return PurificationResult(
        rho,
        method,
        converged,
        termination_reason,
        iterations,
        trace_value,
        abs(trace_value - Ne),
        idempotency_residual,
        hermiticity_residual,
        maxlinkdim(rho),
        spectral_bounds,
        spectral_bounds_validation,
        params.itensors_tol,
        params.itensors_maxdim,
    )
end
