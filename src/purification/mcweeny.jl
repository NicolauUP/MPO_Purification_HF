
"""
    construct_rho_0(sys, H, ϵ, maxχ, H_max, H_min, Ne)

Build the initial density matrix guess by linearly mapping the
eigenvalues of H into [0,1] with the correct electron count Ne.
"""
function construct_rho_0(sys::System, params::ModelParameters ,H_max::Float64, H_min::Float64)
    N = 2^length(sys.sites)
    Ne = round(Int, N * params.density)
    Id = Identity_MPO(sys.sites)   # built internally!
    #= 
    TODO:
    - Add the mean field structure !
    =#
    H = +(sys.H0, sys.VH, sys.VF; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

    μ = real(tr(H) / N) #Technically it should then sum the mean field Hamiltonian. 
    λ = minimum((Ne / (H_max - μ), (N - Ne) / (μ - H_min)))
    coeff_I = (Ne + λ * μ) / N
    coeff_H = -(λ / N)
    return +(coeff_I * Id, coeff_H * H; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
end


"""
    perform_purification(ρ0; maxχ, ϵ, max_steps, verbose)

Adaptive purification scheme. Starts with a trace-correcting linear update
and switches to McWeeny (3P² - 2P³) when idempotency error is small enough.

Returns the purified density matrix ρ, or a partially purified ρ with a
warning if convergence fails.
"""
function perform_purification(ρ0::MPO, params::ModelParameters;verbose::Int=1)
    
    
    use_mcweeny = false
    ρ_old = nothing
    T2_old = 0.0
    idem_error = Inf
    mpo_rel_change = Inf

    for i in 1:params.purification_steps
        if verbose > 0
            println("--- Step $i ---")
        end

        P2 = apply(ρ0, ρ0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        truncate!(P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

        T1 = real(tr(ρ0))
        T2 = real(tr(P2))
        denom = T1 - T2
        idem_error = denom / T1

        if ρ_old !== nothing
            cross_term = real(inner(ρ_old, ρ0))
            diff_norm_sq = max(0.0, T2 + T2_old - 2.0 * cross_term)
            mpo_rel_change = sqrt(diff_norm_sq) / sqrt(T2)
        else
            mpo_rel_change = Inf
        end

        if verbose > 0
            println("  Trace (Ne)           : $T1")
            println("  Rel Idempotency Error (%): $(idem_error * 100)")
            println("  Rel Change in MPO (%)   : $(mpo_rel_change * 100)")
        end

        # Convergence check
        if idem_error < 0.1/100 && mpo_rel_change < 0.1/100
            if verbose > 0 println("\nConverged!\n") end
            print("Final Trace (Ne): $T1, Idempotency Error: $(idem_error * 100) %\n")
            return ρ0
        end

        # Safe break — MPO stuck but not idempotent
        if mpo_rel_change < 1e-8 && idem_error > 0.1/100 # Idempotency error larger than 0.1% but MPO is no longer changing!
            @warn "Purification stuck: MPO is no longer changing (mpo_rel_change < 1e-8) " *
                  "but idempotency error is still large (idem_error = $(idem_error*100) %). " *
                  "Consider increasing maxχ (current: $maxχ)."
            break
        end

        ρ_old = copy(ρ0)
        T2_old = T2

        P3 = apply(ρ0, P2; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        truncate!(P3; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

        if verbose > 0
            println(" χ_1: ", maxlinkdim(ρ0),
                    " χ_2: ", maxlinkdim(P2),
                    " χ_3: ", maxlinkdim(P3))
        end

        T3 = real(tr(P3))

        if abs(denom) > 1e-8
            cn = (T2 - T3) / denom
        else
            use_mcweeny = true
        end

        if !use_mcweeny && (i > 5 && abs(cn - 0.5) < 0.2)
            use_mcweeny = true
            if verbose > 0 println("  Switching to McWeeny purification\n") end
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

        truncate!(ρ0; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end

    @warn "Purification did not converge after $(params.purification_steps) steps. " *
          "Final idempotency error: $idem_error. " *
          "Consider increasing max_steps or maxχ (current: $(params.itensors_maxdim))."
    return ρ0
end
