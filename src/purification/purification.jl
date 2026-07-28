"""
    perform_purification(rho0, params; method=:sp2, ...)

Select a purification backend. The legacy selector `:adaptive_pm_mcweeny` is
accepted as an alias for `:palser_manolopoulos`.
"""
function perform_purification(
    rho0::MPO,
    params::AbstractModelParameters;
    method::Symbol=:sp2,
    verbose::Int=1,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
    sp2_idempotency_tolerance::Real=1e-3,
    sp2_trace_tolerance::Union{Nothing,Real}=nothing,
    chemical_potential::Union{Nothing,Real}=nothing,
    mcweeny_form::Symbol=:standard,
    mcweeny_identity_mpo::Union{Nothing,MPO}=nothing,
    mcweeny_trace_target::Union{Nothing,Real}=nothing,
    mcweeny_trace_tolerance::Union{Nothing,Real}=nothing,
    gc_policy::Symbol=:automatic,
    gc_period::Integer=10,
    gc_threshold_bytes::Integer=1 << 30,
    device_cleanup::Function=() -> nothing,
)
    if method == :sp2
        return perform_purification_sp2(
            rho0, params;
            verbose=verbose,
            io=io,
            overwrite_progress=overwrite_progress,
            idempotency_tolerance=sp2_idempotency_tolerance,
            trace_tolerance=sp2_trace_tolerance,
            spectral_bounds=spectral_bounds,
            spectral_bounds_validation=spectral_bounds_validation,
        )
    elseif method == :mcweeny_mu
        return perform_purification_mcweeny_mu(rho0, params;
            verbose=verbose, spectral_bounds=spectral_bounds,
            spectral_bounds_validation=spectral_bounds_validation,
            chemical_potential=chemical_potential,
            form=mcweeny_form,
            identity_mpo=mcweeny_identity_mpo,
            trace_target=mcweeny_trace_target,
            trace_tolerance=mcweeny_trace_tolerance,
        )
    elseif method == :palser_manolopoulos || method == :adaptive_pm_mcweeny
        return perform_purification_palser_manolopoulos(
            rho0, params;
            verbose=verbose,
            io=io,
            overwrite_progress=overwrite_progress,
            spectral_bounds=spectral_bounds,
            spectral_bounds_validation=spectral_bounds_validation,
            gc_policy=gc_policy,
            gc_period=gc_period,
            gc_threshold_bytes=gc_threshold_bytes,
            device_cleanup=device_cleanup,
        )
    end
    throw(ArgumentError(
        "unknown purification method $method; supported methods are :palser_manolopoulos, :mcweeny_mu, and :sp2",
    ))
end
