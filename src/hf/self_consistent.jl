

function run_scf!(sys::System, params::SCFParams, H_min::Float64, H_max::Float64; verbose::Int=:nothing) 
    if verbose == :all
        println("Starting SCF iterations with parameters:")
        println("  Max Iterations: $(params.max_iterations)")
        println("  Convergence Tolerance: $(params.scf_tol)")
        println("  Mixing Parameter: $(params.scf_mixing)")
    end


    converged = false
    rel_change = Inf
    for iter in 1:params.max_iterations
        if verbose == :all
            println("SCF Iteration $iter")

        end
        # Step 1: Obtain density matrix!
        ρ0 = construct_ρ0(sys, params, H_min, H_max)
        ρ_purified = perform_purification(ρ0, params; verbose=verbose == :all ? 1 : 0)
        sys.ρ = ρ_purified
        # Step 2: Extract Hartree potential
        vh_mpo = extract_hartree_mpo_1d(sys)
        

        # Step 3: Mix with previous potential
        if iter == 1
            sys.VH = vh_mpo
        else
            sys.VH = +(params.scf_mixing * vh_mpo, (1 - params.scf_mixing) * sys.VH; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        end 

        # Step 4: Check convergence        
        if iter > 1
            vh_diff = +(vh_mpo, -sys.VH; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
            diff_norm = sqrt(real(inner(vh_diff, vh_diff)))
            vh_norm = sqrt(real(inner(sys.VH, sys.VH)))
            rel_change = diff_norm / vh_norm    
        end
        if verbose == :all
            println("  Relative Change in VH: $(rel_change * 100) %")
        end
        if iter > 1 && rel_change < params.scf_tol
            converged = true
            if verbose == :all
                println("\nSCF Converged in $iter iterations with relative change $(rel_change * 100) %\n")
            end
            break
        end
    end
    if !converged && verbose != :nothing
        println("\nSCF did not converge within $(params.max_iterations) iterations. Final relative change: $(rel_change * 100) %\n")
    end
    return converged    

end