

function run_scf!(sys::System, params::SCFParams, H_min::Float64, H_max::Float64; verbose::Int=:nothing) 
    if verbose == :all
        println("Starting SCF iterations with parameters:")
        println("  Max Iterations: $(params.max_iterations)")
        println("  Convergence Tolerance: $(params.convergence_tol)")
        println("  Mixing Parameter: $(params.mixing_parameter)")
    end

    converged = false
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
        

        
    end
end