using CUDA
"""
    construct_rho_0(sys, H, ϵ, maxχ, H_max, H_min, Ne)

Build the initial density matrix guess by linearly mapping the
eigenvalues of H into [0,1] with the correct electron count Ne.
"""
function construct_rho_0(sys::System, params::AbstractModelParameters ,H_min::Float64, H_max::Float64;
    to_gpu=identity)
    
    N = 2^length(sys.sites)
    Ne = round(Int, N * params.density)
    Id_cpu = Identity_MPO(sys.sites)   # built internally!
    Id = to_gpu(Id_cpu) 

    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim) #Computes the MF-Hamiltonian

    μ = real(tr(H) / N) 
    @assert H_max > μ > H_min "μ=$μ outside [H_min=$H_min, H_max=$H_max]" 

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

Returns the purified density matrix ρ, or a partially purified ρ with a
warning if convergence fails.
"""
function perform_purification(ρ0::MPO, params::AbstractModelParameters;verbose::Int=1)

    N = 2^params.L
    println("N = $N, density = $(params.density), Ne = $(round(Int, N * params.density))")
    Ne = round(Int, N * params.density)
 
    cn = nothing
    use_mcweeny = false
    T2_old = 0.0
    idem_error = Inf


    for i in 1:params.purification_steps
        if verbose > 0
            println("--- Step $i ---")
        end

        P2 = apply(ρ0, ρ0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)


        T1 = real(tr(ρ0))
        T2 = real(tr(P2))
        denom = T1 - T2 

        idem_error = denom / T1 

        if verbose > 0
            println("  Trace (Ne)           : $T1")
            println("  Drift in Trace (%)  : $((T1 - Ne) / Ne * 100)")
            println("  Rel Idempotency Error (%): $(idem_error * 100)")
        end
      



        P3 = apply(ρ0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

        if verbose > 0
            println(" χ_1: ", maxlinkdim(ρ0),"\n",
                    " χ_2: ", maxlinkdim(P2),"\n",
                    " χ_3: ", maxlinkdim(P3))
        end

        T3 = real(tr(P3))
        

       if abs(T1 - Ne) / Ne > 0.1/100
            @warn "Trace has drifted: T1=$T1, Ne=$Ne. Stopping purification."
            return ρ0
       end
       
       if abs(denom) < 1e-3 #Hard cutoff to avoid numerical issues - test!
            use_mcweeny = true
       else
            cn = (T2 - T3) / denom  # normal PM step
            if abs(cn - 0.5)/0.5 * 100 < 0.1 # 0.1% from 0.5.
                use_mcweeny = true  
            end
        end

        if abs(denom)/T1*100 < 1e-1  #0.1# Error in the idempotency is small enough,  Stopping!
            @info "T1 Similar to T2, converged!"
            return ρ0
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

    @warn "Purification did not converge after $(params.purification_steps) steps. " *
          "Final idempotency error: $idem_error. " *
          "Consider increasing max_steps or maxχ (current: $(params.itensors_maxdim))."
    return ρ0
end


function perform_purification_grandcanonical(sys::System, params::AbstractModelParameters,H_min::Float64,H_max::Float64;verbose::Int=1, to_gpu=identity)

    N = 2^params.L
    Ne = round(Int, N * params.density)

    if verbose > 0
        println("N = $N, density = $(params.density), Ne = $Ne")
    end

    idempotency_tol = 1e-3 # 0.1% idempotency error tolerance
    #Here, we will implement a bissection search of the chemical potential, to avoid precise trace calculations, which fail with bound dimension limits!


    mu_low = H_min
    mu_high = H_max
    rho_new = nothing #Outer scope 

    while (mu_high - mu_low) > 1e-7
        mu = (mu_high + mu_low) / 2.0
        
        if verbose > 0
            println("Current μ: $mu, μ_low: $mu_low, μ_high: $mu_high")
        end

        rho_0 = construct_rho_0_mcweeny(sys, params, mu, H_min, H_max; to_gpu=to_gpu)
        
        abort_flag = false

        for i in 1:params.purification_steps
            if verbose > 0 
                println("     --- Step $i ---")
            end
            T0 = real(tr(rho_0))        
            println("     Trace (Ne)           : $T0")

            P2 = apply(rho_0, rho_0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            T2 = real(tr(P2))

            idem_error = abs(T2 - T0) / max(T0,1.0)  # Avoid division by zero
            if idem_error < idempotency_tol

                if verbose > 0 
                    println("Converged at step $i. Trace: $T0, Idempotency Error: $idem_error")
                end
                return rho_0
            end




            P3 = apply(rho_0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            rho_new = +(3.0 * P2, -2.0 * P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim) #Simple Mcweeny Update!
            println(" Bond: ", maxlinkdim(rho_new))
            T1 = real(tr(rho_new))

            #======== Verify if we are above or below the target Ne =#
            if T1 - Ne > 0.5 && (T1 - T0) > 0  #We are above the target, by 0.5 particles and the trace is increasing!
                mu_high = mu
                abort_flag = true
                break
            elseif T1 - Ne < -0.5 && (T1 - T0) < 0  #We are below the target, by 0.5 particles and the trace is decreasing!
                mu_low = mu
                abort_flag = true
                break
            end

            rho_0 = rho_new
        end

        if !abort_flag
            @warn "Inner loop maxed out without aborting. Forcing bissection update."
            if real(tr(rho_new)) > Ne
                mu_high = mu
            else
                mu_low = mu
            end
        end
    end
    @warn "Bissection resolution limit reached."
    return rho_new
end

            


