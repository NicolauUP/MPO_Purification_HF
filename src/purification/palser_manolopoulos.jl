"""
    perform_purification_palser_manolopoulos(ρ0, params; verbose)

Adaptive purification scheme. Starts with a trace-correcting linear update
and switches to McWeeny (3P² - 2P³) when idempotency error is small enough.

Returns a [`PurificationResult`](@ref) containing the density matrix and
termination diagnostics. The current adaptive PM/McWeeny stopping criterion
is preserved unchanged.
"""
function perform_purification_palser_manolopoulos(
    ρ0::MPO,
    params::AbstractModelParameters;
    verbose::Int=1,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
    gc_policy::Symbol=:automatic,
    gc_period::Integer=10,
    gc_threshold_bytes::Integer=1 << 30,
    device_cleanup::Function=() -> nothing,
)

    _validate_gc_policy(gc_policy, gc_period, gc_threshold_bytes)

    if !isnothing(spectral_bounds)
        spectral_bounds = validate_spectral_bounds(spectral_bounds...)
    end

    N = 2^params.L
    Ne = round(Int, N * params.density)
    verbose > 0 && println(io, "Purifying N=$N, density=$(params.density), Ne=$Ne")
 
    cn = nothing
    use_mcweeny = false
    T2_old = 0.0
    idem_error = Inf
    idempotency_tolerance = 1e-3
    max_bond_dimension = maxlinkdim(ρ0)
    bond_dimension_sum = 0
    bond_dimension_samples = 0
    bond_dimensions = NTuple{3,Int}[]


    for i in 1:params.purification_steps
        P2 = apply(ρ0, ρ0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)


        T1 = real(tr(ρ0))
        T2 = real(tr(P2))
        denom = T1 - T2 

        idem_error = idempotency_residual(ρ0, P2)

        P3 = apply(ρ0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        max_bond_dimension = max(max_bond_dimension, maxlinkdim(ρ0), maxlinkdim(P2), maxlinkdim(P3))
        bond_dimension_sum += maxlinkdim(ρ0) + maxlinkdim(P2) + maxlinkdim(P3)
        bond_dimension_samples += 3
        push!(bond_dimensions, (maxlinkdim(ρ0), maxlinkdim(P2), maxlinkdim(P3)))
        work = PurificationWorkStats(
            i, i, max_bond_dimension, bond_dimension_sum / bond_dimension_samples,
            copy(bond_dimensions),
        )

        if verbose > 0
            details = @sprintf(
                "Tr=%.12g | drift=%+.3e%% | idem=%.3e%% | χ=(%d,%d,%d)",
                T1,
                (T1 - Ne) / Ne * 100,
                idem_error * 100,
                maxlinkdim(ρ0),
                maxlinkdim(P2),
                maxlinkdim(P3),
            )
            print_iteration_progress(
                io, "Purification", i, params.purification_steps, details;
                overwrite=overwrite_progress,
            )
        end

        T3 = real(tr(P3))
        

       if abs(T1 - Ne) / Ne > 0.1/100
            verbose > 0 && finish_iteration_progress(io, overwrite_progress)
            @warn "Trace has drifted: T1=$T1, Ne=$Ne. Stopping purification."
            return purification_result(
                ρ0, params;
                method=:palser_manolopoulos,
                converged=false,
                termination_reason=:trace_drift,
                iterations=i,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
                work=work,
            )
       end
       
       if abs(denom) < 1e-3 #Hard cutoff to avoid numerical issues - test!
            use_mcweeny = true
       else
            cn = (T2 - T3) / denom  # normal PM step
            if abs(cn - 0.5)/0.5 * 100 < 0.1 # 0.1% from 0.5.
                use_mcweeny = true  
            end
        end

        if idem_error < idempotency_tolerance
            verbose > 0 && finish_iteration_progress(io, overwrite_progress)
            verbose > 0 && println(io, "Purification converged at step $i.")
            return purification_result(
                ρ0, params;
                method=:palser_manolopoulos,
                converged=true,
                termination_reason=:idempotency_threshold,
                iterations=i,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
                work=work,
            )
        end
        if use_mcweeny
            ρ0 = +(3.0 * P2, -2.0 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        else
            if cn < 0.5
                inv_fac = 1.0 / (1.0 - cn)
                c1 = (1 - 2cn) * inv_fac
                c2 = (1 + cn) * inv_fac
                c3 = -1.0 * inv_fac
                ρ0 = +(c1 * ρ0, c2 * P2, c3 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            else
                inv_fac = 1.0 / cn
                c2 = (1.0 + cn) * inv_fac
                c3 = -1.0 * inv_fac
                ρ0 = +(c2 * P2, c3 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            end
        end


 # --- Memory Management ---
        P2 = nothing
        P3 = nothing
        
        maybe_collect_garbage!(i;
            gc_policy=gc_policy,
            gc_period=gc_period,
            gc_threshold_bytes=gc_threshold_bytes,
        )
        device_cleanup()

    end

    verbose > 0 && finish_iteration_progress(io, overwrite_progress)
    @warn "Purification did not converge after $(params.purification_steps) steps. " *
          "Final normalized idempotency residual: $idem_error. " *
          "Consider increasing max_steps or maxχ (current: $(params.itensors_maxdim))."
    return purification_result(
        ρ0, params;
        method=:palser_manolopoulos,
        converged=false,
        termination_reason=:max_iterations,
        iterations=params.purification_steps,
        spectral_bounds=spectral_bounds,
        spectral_bounds_validation=spectral_bounds_validation,
        work=PurificationWorkStats(
            params.purification_steps, params.purification_steps, max_bond_dimension,
            bond_dimension_sum / max(1, bond_dimension_samples),
            copy(bond_dimensions),
        ),
    )
end
