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


#= 
TODO:
    1) Maybe I can test that version that only needs rho and rho^2 instead of rho^3? It requires more iterations but maybe it is more stable?   
=#
