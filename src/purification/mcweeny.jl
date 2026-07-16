using CUDA
"""
    construct_rho_0(sys, params, H_min, H_max; verify_spectral_bounds=false)

Build the initial density matrix guess by linearly mapping the
effective Hamiltonian into [0,1] with the correct electron count. Bounds are
user supplied. Setting `verify_spectral_bounds=true` performs an exact,
small-system CPU validation before scaling.
"""
function construct_rho_0(
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
    Id_cpu = Identity_MPO(sys.sites)   # built internally!
    Id = to_gpu(Id_cpu) 

    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim) #Computes the MF-Hamiltonian

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
    
    ρ_temp = +(coeff_I * Id, coeff_H * H; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)


    Id = nothing
    Id_cpu = nothing 
    return ρ_temp
end

function construct_rho_0_mcweeny(sys::System, params::AbstractModelParameters, mu, H_min::Float64, H_max::Float64; to_gpu=identity)
    N = 2^length(sys.sites)
    Ne = round(Int, N * params.density)
    Id_cpu = Identity_MPO(sys.sites)   # built internally!
    Id = to_gpu(Id_cpu) 

    c = 0.5 / max(mu - H_min, H_max - mu)
    coeff_I = 0.5 + c * mu
    coeff_H = -c
    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim) #Computes the MF-Hamiltonian
    rho_0 = +(coeff_I * Id, coeff_H * H; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return rho_0
end

"""
    perform_purification(ρ0; maxχ, ϵ, max_steps, verbose)

Adaptive purification scheme. Starts with a trace-correcting linear update
and switches to McWeeny (3P² - 2P³) when idempotency error is small enough.

Returns a [`PurificationResult`](@ref) containing the density matrix and
termination diagnostics. The current adaptive PM/McWeeny stopping criterion
is preserved unchanged.
"""
function perform_purification(
    ρ0::MPO,
    params::AbstractModelParameters;
    verbose::Int=1,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    spectral_bounds::Union{Nothing,Tuple{Float64,Float64}}=nothing,
    spectral_bounds_validation::Symbol=:not_provided,
)

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


    for i in 1:params.purification_steps
        P2 = apply(ρ0, ρ0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)


        T1 = real(tr(ρ0))
        T2 = real(tr(P2))
        denom = T1 - T2 

        idem_error = idempotency_residual(ρ0, P2)

        P3 = apply(ρ0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

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
                method=:adaptive_pm_mcweeny,
                converged=false,
                termination_reason=:trace_drift,
                iterations=i,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
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
                method=:adaptive_pm_mcweeny,
                converged=true,
                termination_reason=:idempotency_threshold,
                iterations=i,
                spectral_bounds=spectral_bounds,
                spectral_bounds_validation=spectral_bounds_validation,
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
        
        # Force garbage collection
        GC.gc()
        if CUDA.functional()
            CUDA.reclaim()  # Reclaim GPU memory if using CUDA
        end

    end

    verbose > 0 && finish_iteration_progress(io, overwrite_progress)
    @warn "Purification did not converge after $(params.purification_steps) steps. " *
          "Final normalized idempotency residual: $idem_error. " *
          "Consider increasing max_steps or maxχ (current: $(params.itensors_maxdim))."
    return purification_result(
        ρ0, params;
        method=:adaptive_pm_mcweeny,
        converged=false,
        termination_reason=:max_iterations,
        iterations=params.purification_steps,
        spectral_bounds=spectral_bounds,
        spectral_bounds_validation=spectral_bounds_validation,
    )
end
function perform_purification_grandcanonical(sys::System, params::AbstractModelParameters, H_min::Float64, H_max::Float64; verbose::Int=1, to_gpu=identity)

    N = 2^params.L
    Ne = round(Int, N * params.density)

    if verbose > 0
        println("N = $N, density = $(params.density), Ne = $Ne")
    end

    idempotency_tol = 1e-3 
    mu_low = -1.0
    mu_high = 1.0
    rho_new = nothing 

    while (mu_high - mu_low) > 1e-4
        mu = (mu_high + mu_low) / 2.0
        
        if verbose > 0
            println("Current μ: $mu, μ_low: $mu_low, μ_high: $mu_high")
        end

        rho_0 = construct_rho_0_mcweeny(sys, params, mu, H_min, H_max; to_gpu=to_gpu)
        
        for i in 1:params.purification_steps
            if verbose > 0 
                println("     --- Step $i ---")
            end
            
            T0 = real(tr(rho_0))        
            
            if verbose > 0
                println("     Trace (Ne)           : $T0")
            end

            P2 = apply(rho_0, rho_0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            T2 = real(tr(P2))

            P3 = apply(rho_0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            rho_new = +(3.0 * P2, -2.0 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim) 
            
            if verbose > 0
                println(" Bond: ", maxlinkdim(rho_new))
            end

            T1 = real(tr(rho_new))
            idem_error = abs(T2 - T0) / max(T0, 1.0)  
            
            if idem_error < idempotency_tol
                if verbose > 0 
                    println("Converged at step $i. Trace: $T0, Idempotency Error: $idem_error")
                end
                break
            end
            rho_0 = rho_new
        end

        T_final = real(tr(rho_new))

        if abs(T_final - Ne) < 0.5
            if verbose > 0 
                println("Final trace $T_final is within tolerance of Ne=$Ne. Converged.")
            end
            return rho_new
        elseif T_final > Ne
            mu_high = mu
        else
            mu_low = mu
        end
    end # Added missing end statement

    @warn "Bisection resolution limit reached."
    return rho_new
end
