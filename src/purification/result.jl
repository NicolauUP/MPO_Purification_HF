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
    truncation_cutoff::Float64
    maxdim::Int
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
        params.itensors_tol,
        params.itensors_maxdim,
    )
end
