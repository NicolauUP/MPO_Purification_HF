

function run_scf!(sys::System, H_min::Float64, H_max::Float64; 
    verbose::Symbol=:nothing,
    to_gpu=identity,
    to_cpu=identity,
    cleanup= () -> nothing)

    if verbose == :all
        println("Starting SCF iterations with parameters:")
        println("  Max Iterations: $(sys.params.scf_max_iterations)")
        println("  Convergence Tolerance: $(sys.params.scf_tol)")
        println("  Mixing Parameter: $(sys.params.scf_mixing)")
    elseif verbose == :Density
        println("Starting SCF iterations with parameters:")
        println("  Max Iterations: $(sys.params.scf_max_iterations)")
        println("  Convergence Tolerance: $(sys.params.scf_tol)")
        println("  Mixing Parameter: $(sys.params.scf_mixing)")
    end

    sys.H0 = to_gpu(sys.H0)
    sys.VH = to_gpu(sys.VH)
    sys.VF = to_gpu(sys.VF) #

   
    params = sys.params

    converged = false
    rel_change = Inf
    for iter in 1:params.scf_max_iterations
        if verbose != :nothing
            println("\n ----------- SCF Iteration $iter")

        end
        # Step 1: Obtain density matrix!
        ρ0_device = construct_rho_0(sys, params, H_min, H_max; to_gpu=to_gpu)
        if verbose == :Density
            T1 = real(tr(ρ0_device))
            println("  Trace (Ne) of initial ρ0: $T1")
        end

        ρ_purified_device = perform_purification(ρ0_device, params; verbose=1)

        if verbose == :Density
            T1_purified = real(tr(ρ_purified_device))
            println("  Trace (Ne) of purified ρ: $T1_purified")
        end
        sys.ρ = to_cpu(ρ_purified_device) #Ok updates rho!
        #=
        Verify how this to_cpu should be done! 
        There are two options:
        1. Move the purified density matrix back to CPU and then extract the Hartree potential
        2. Extract the Hartree potential directly on the GPU and only move the resulting MPO back to CPU. 

        Verify which one is more efficient. I believe 1 is ok, because TCI is scalar! 
        TODO: Write to_cpu funcion! 
        =#

        # Step 2: Extract Hartree potential
        vh_mpo_cpu = extract_hartree_mpo_1d(sys)
        vf_mpo_cpu = extract_fock_mpo_1d(sys) 

        vh_mpo = to_gpu(vh_mpo_cpu)
        vf_mpo = to_gpu(vf_mpo_cpu)



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

        if iter == 1
            sys.VH = vh_mpo
            sys.VF = vf_mpo
        else
            sys.VH = +(sys.params.scf_mixing * vh_mpo, (1 - params.scf_mixing) * sys.VH; cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim)
            sys.VF = +(sys.params.scf_mixing * vf_mpo, (1 - params.scf_mixing) * sys.VF; cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim)
        end 

        if iter > 1 && rel_change * 100 < params.scf_tol
            converged = true
            if verbose == :all
                println("\nSCF Converged in $iter iterations with relative change $(rel_change * 100) %\n")
            end
            sys.VH = vh_mpo # Final update to ensure we return the most accurate VH
            sys.VF = vf_mpo # Final update to ensure we return the most accurate VF
            break
        end
        # Optional cleanup to free GPU memory after each iteration
        cleanup() # Default to a no-op, but can be set to a function that frees GPU memory if needed!
        GC.gc(true) # Force garbage collection to free memory
    end
    if !converged && verbose != :nothing
        println("\nSCF did not converge within $(params.scf_max_iterations) iterations. Final relative change: $(rel_change * 100) %\n")
    end
    return converged    

end