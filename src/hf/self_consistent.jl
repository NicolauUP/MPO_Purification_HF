function _mpo_relative_change(candidate::MPO, reference::MPO, params::AbstractModelParameters)
    difference = +(candidate, -reference; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return safe_relative_change(
        sqrt(max(0.0, real(inner(difference, difference)))),
        sqrt(max(0.0, real(inner(reference, reference)))),
    )
end

function _scf_commutator_residual(H::MPO, rho::MPO, params::AbstractModelParameters)
    Hrho = apply(H, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rhoH = apply(rho, H; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return _mpo_relative_change(Hrho, rhoH, params)
end

function run_scf!(sys::System, H_min::Float64, H_max::Float64; 
    verbose::Symbol=:nothing,
    allow_unconverged_purification::Bool=false,
    verify_spectral_bounds::Bool=false,
    spectral_safety_margin::Float64=0.0,
    purification_method::Symbol=:sp2,
    chemical_potential::Union{Nothing,Real}=nothing,
    to_gpu=identity,
    to_cpu=identity,
    cleanup= () -> nothing)

    sys.params isa ParametersSquare && throw(ArgumentError(
        "square-lattice SCF is unsupported: square Hartree/Fock field extraction is not implemented",
    ))


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
    vh_residual = Inf
    vf_residual = Inf
    rho_residual = Inf
    commutator_residual = Inf
    for iter in 1:params.scf_max_iterations
        if verbose != :nothing
            println("\n ----------- SCF Iteration $iter")

        end
        # Step 1: Obtain density matrix!
        ρ0_device = construct_rho_0(
            sys, params, H_min, H_max;
            to_gpu=to_gpu,
            verify_spectral_bounds=verify_spectral_bounds,
            safety_margin=spectral_safety_margin,
            method=purification_method,
            chemical_potential=chemical_potential,
        )
        if verbose == :Density
            T1 = real(tr(ρ0_device))
            println("  Trace (Ne) of initial ρ0: $T1")
        end

        purification_verbose = verbose == :nothing ? 0 : 1
        purification = perform_purification(
            ρ0_device, params;
            verbose=purification_verbose,
            spectral_bounds=(H_min, H_max),
            spectral_bounds_validation=(
                verify_spectral_bounds ? :exact_small_system : :supplied_unverified
            ),
            method=purification_method,
            chemical_potential=chemical_potential,
        )
        if !purification.converged && !allow_unconverged_purification
            verbose != :nothing && println(
                "SCF stopped: purification ended with $(purification.termination_reason).",
            )
            return false
        end
        ρ_purified_device = purification.rho

        if verbose == :Density
            T1_purified = real(tr(ρ_purified_device))
            println("  Trace (Ne) of purified ρ: $T1_purified")
        end
        rho_previous = sys.ρ
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
        H_effective = +(sys.H0, sys.VH, sys.VF;
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
        )
        commutator_residual = _scf_commutator_residual(H_effective, ρ_purified_device, params)
        if iter > 1
            vh_residual = _mpo_relative_change(vh_mpo, sys.VH, params)
            vf_residual = _mpo_relative_change(vf_mpo, sys.VF, params)
            rho_residual = _mpo_relative_change(sys.ρ, rho_previous, params)
        end
        if verbose == :all
            println("  Residuals (%): VH=$(vh_residual * 100), VF=$(vf_residual * 100), ρ=$(rho_residual * 100), [H,ρ]=$(commutator_residual * 100)")
        end

        if iter > 1 && maximum((vh_residual, vf_residual, rho_residual, commutator_residual)) * 100 < params.scf_tol
            converged = true
            if verbose == :all
                println("\nSCF converged in $iter iterations.\n")
            end
            break
        end

        if iter == 1
            sys.VH = vh_mpo
            sys.VF = vf_mpo
        else
            sys.VH = +(sys.params.scf_mixing * vh_mpo, (1 - params.scf_mixing) * sys.VH; cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim)
            sys.VF = +(sys.params.scf_mixing * vf_mpo, (1 - params.scf_mixing) * sys.VF; cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim)
        end
        # Optional cleanup to free GPU memory after each iteration
        cleanup() # Default to a no-op, but can be set to a function that frees GPU memory if needed!
        GC.gc(true) # Force garbage collection to free memory
    end
    if !converged && verbose != :nothing
        println("\nSCF did not converge within $(params.scf_max_iterations) iterations. Final maximum residual: $(maximum((vh_residual, vf_residual, rho_residual, commutator_residual)) * 100) %\n")
    end
    return converged    

end
