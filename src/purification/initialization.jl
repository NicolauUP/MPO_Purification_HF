"""
    construct_rho_0(sys, params, H_min, H_max;
                    method=:palser_manolopoulos, verify_spectral_bounds=false)

Build an initial density-matrix MPO from the effective Hamiltonian.
`method=:palser_manolopoulos` uses the trace-correcting linear map, while
`method=:sp2` uses the canonical SP2 spectral map. Bounds are user supplied.
Setting `verify_spectral_bounds=true` performs an exact, small-system CPU
validation before scaling.
"""
function construct_rho_0(
    sys::System,
    params::AbstractModelParameters,
    H_min::Float64,
    H_max::Float64;
    to_gpu=identity,
    verify_spectral_bounds::Bool=false,
    safety_margin::Float64=0.0,
    method::Symbol=:palser_manolopoulos,
    chemical_potential::Union{Nothing,Real}=nothing,
)
    if method == :palser_manolopoulos || method == :adaptive_pm_mcweeny
        return _construct_rho_0_adaptive(
            sys, params, H_min, H_max;
            to_gpu=to_gpu,
            verify_spectral_bounds=verify_spectral_bounds,
            safety_margin=safety_margin,
        )
    elseif method == :sp2
        return _construct_rho_0_sp2(
            sys, params, H_min, H_max;
            to_gpu=to_gpu,
            verify_spectral_bounds=verify_spectral_bounds,
            safety_margin=safety_margin,
        )
    elseif method == :mcweeny_mu
        isnothing(chemical_potential) && throw(ArgumentError("method=:mcweeny_mu requires chemical_potential"))
        return construct_rho_0_mcweeny(sys, params, chemical_potential, H_min, H_max; to_gpu=to_gpu)
    end
    throw(ArgumentError(
        "unknown purification method $method; supported methods are :palser_manolopoulos, :mcweeny_mu, and :sp2",
    ))
end

function _construct_rho_0_adaptive(
    sys::System,
    params::AbstractModelParameters,
    H_min::Float64,
    H_max::Float64;
    to_gpu=identity,
    verify_spectral_bounds::Bool=false,
    safety_margin::Float64=0.0,
)
    N = 2^length(sys.sites)
    Ne = round(Int, N * params.density)
    Id = to_gpu(Identity_MPO(sys.sites))
    H = +(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )

    validate_spectral_bounds(H_min, H_max; safety_margin=safety_margin)
    verify_spectral_bounds && verify_spectral_bounds_exact(
        sys, H, H_min, H_max; safety_margin=safety_margin,
    )
    μ = real(tr(H) / N)
    H_max > μ > H_min || throw(ArgumentError(
        "mean energy μ=$μ lies outside supplied spectral bounds [$H_min, $H_max]",
    ))

    λ = minimum((Ne / (H_max - μ), (N - Ne) / (μ - H_min)))
    coeff_I = (Ne + λ * μ) / N
    coeff_H = -(λ / N)
    return +(coeff_I * Id, coeff_H * H;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

function _construct_rho_0_sp2(
    sys::System,
    params::AbstractModelParameters,
    H_min::Float64,
    H_max::Float64;
    to_gpu=identity,
    verify_spectral_bounds::Bool=false,
    safety_margin::Float64=0.0,
)
    Id = to_gpu(Identity_MPO(sys.sites))
    H = +(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )

    validate_spectral_bounds(H_min, H_max; safety_margin=safety_margin)
    verify_spectral_bounds && verify_spectral_bounds_exact(
        sys, H, H_min, H_max; safety_margin=safety_margin,
    )
    width = H_max - H_min
    return +(
        (H_max / width) * Id,
        (-1.0 / width) * H;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
end
