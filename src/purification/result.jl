"""
    PurificationResult

Diagnostics returned by a density-matrix purification backend. `spectral_bounds`
is `nothing` when the backend was not given explicit bounds.
"""
struct PurificationWorkStats
    squares::Int
    cubes::Int
    max_bond_dimension::Int
    mean_bond_dimension::Float64
    bond_dimensions::Vector{NTuple{3,Int}}
end

struct PurificationResult{R}
    rho::R
    method::Symbol
    converged::Bool
    termination_reason::Symbol
    iterations::Int
    selected_iteration::Int
    trace::Float64
    target_particles::Union{Nothing,Float64}
    trace_error::Union{Nothing,Float64}
    idempotency_residual::Float64
    hermiticity_residual::Float64
    final_bond_dimension::Int
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}
    spectral_bounds_validation::Symbol
    truncation_cutoff::Float64
    maxdim::Int
    work::PurificationWorkStats
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
    square_scf_spectral_bounds(params; potential_bounds=nothing,
                               hopping_abs_bounds=nothing, margin=0.0)

Return a conservative interval for every effective Hamiltonian of the open
square nearest-neighbour Hartree--Fock model, provided the supplied
coefficient bounds are valid and the density matrix remains a contraction.

The open square has coordination at most four. The implemented Hartree and
real-exchange fields therefore have row-sum bounds `4abs(U)` each. Together
with the hopping row sum, the returned interval is

```math
[W_{min} - 2(t_x^* + t_y^*) - 8|U| - m,
 W_{max} + 2(t_x^* + t_y^*) + 8|U| + m].
```

Numeric hopping bounds are inferred from `params.t`. Functional hopping needs
`hopping_abs_bounds=(max_abs_tx, max_abs_ty)`. A functional external potential
needs its known global range as `potential_bounds=(W_min, W_max)`. `margin` is
additional nonnegative energy padding for MPO/truncation uncertainty.
"""
function square_scf_spectral_bounds(
    params::ParametersSquare;
    potential_bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    hopping_abs_bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
    margin::Real=0.0,
)
    isfinite(margin) && margin >= 0 || throw(ArgumentError(
        "margin must be finite and nonnegative, got $margin",
    ))
    hopping_bounds = if isnothing(hopping_abs_bounds)
        all(component -> component isa Number, params.t) || throw(ArgumentError(
            "functional square hopping requires hopping_abs_bounds=(abs_tx, abs_ty)",
        ))
        (abs(Float64(params.t[1])), abs(Float64(params.t[2])))
    else
        hopping_abs_bounds
    end
    all(value -> isfinite(value) && value >= 0, hopping_bounds) || throw(ArgumentError(
        "hopping_abs_bounds must be finite and nonnegative, got $hopping_bounds",
    ))
    onsite_bounds = if isnothing(potential_bounds)
        isnothing(params.W) || throw(ArgumentError(
            "a square external potential function requires potential_bounds=(W_min, W_max)",
        ))
        (0.0, 0.0)
    else
        potential_bounds
    end
    W_min, W_max = Float64.(onsite_bounds)
    isfinite(W_min) && isfinite(W_max) && W_min <= W_max || throw(ArgumentError(
        "potential_bounds must be finite and satisfy W_min <= W_max, got $onsite_bounds",
    ))
    hopping_radius = 2 * (Float64(hopping_bounds[1]) + Float64(hopping_bounds[2]))
    interaction_radius = 8 * abs(Float64(params.U))
    return validate_spectral_bounds(
        W_min - hopping_radius - interaction_radius - margin,
        W_max + hopping_radius + interaction_radius + margin,
    )
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

function idempotency_residual(rho::MPO, rho_squared::MPO)
    rho_norm_squared = max(0.0, real(inner(rho, rho)))
    difference_norm_squared = max(
        0.0,
        real(inner(rho_squared, rho_squared)) + rho_norm_squared -
        2 * real(inner(rho_squared, rho)),
    )
    return sqrt(difference_norm_squared) / max(sqrt(rho_norm_squared), sqrt(eps(Float64)))
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
    target_particles::Union{Nothing,Real}=round(Int, 2^params.L * params.density),
    selected_iteration::Integer=iterations,
    work::PurificationWorkStats=PurificationWorkStats(
        0, 0, maxlinkdim(rho), Float64(maxlinkdim(rho)), [(maxlinkdim(rho), 0, 0)],
    ),
)
    1 <= selected_iteration <= iterations || throw(ArgumentError(
        "selected_iteration must satisfy 1 <= selected_iteration <= iterations",
    ))
    rho_squared = apply(rho, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    trace_value = Float64(real(tr(rho)))
    idem_residual = idempotency_residual(rho, rho_squared)
    hermiticity_residual = _relative_mpo_residual(rho, ITensors.dag(rho), params)

    return PurificationResult(
        rho,
        method,
        converged,
        termination_reason,
        iterations,
        Int(selected_iteration),
        trace_value,
        isnothing(target_particles) ? nothing : Float64(target_particles),
        isnothing(target_particles) ? nothing : abs(trace_value - target_particles),
        idem_residual,
        hermiticity_residual,
        maxlinkdim(rho),
        spectral_bounds,
        spectral_bounds_validation,
        params.itensors_tol,
        params.itensors_maxdim,
        work,
    )
end
