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

function _scf_energy_total(sys::System)
    if sys.params isa Parameters1D
        return nearest_neighbor_hf_energy_1d(sys).total
    elseif sys.params isa ParametersSquare
        return nearest_neighbor_hf_energy_square(sys).total
    end
    return nothing
end

function _scf_record_within_tolerance(record::SCFIterationRecord, tolerance::Real)
    residuals = (
        record.vh_residual,
        record.vf_residual,
        record.rho_residual,
        record.commutator_residual,
    )
    return all(isfinite, residuals) && maximum(residuals) <= tolerance
end

function _scf_has_stable_history(
    history::Vector{SCFIterationRecord},
    tolerance::Real;
    required_iterations::Integer=2,
)
    required_iterations > 0 || throw(ArgumentError("required_iterations must be positive"))
    length(history) >= required_iterations || return false
    return all(record -> _scf_record_within_tolerance(record, tolerance),
        @view history[(end - required_iterations + 1):end])
end

function _scf_is_two_cycle(record::SCFIterationRecord, tolerance::Real)
    return isfinite(record.two_cycle_residual) &&
           record.two_cycle_residual <= tolerance &&
           isfinite(record.rho_residual) && record.rho_residual > tolerance
end

function run_scf!(sys::System, H_min::Float64, H_max::Float64; 
    verbose::Symbol=:nothing,
    allow_unconverged_purification::Bool=false,
    verify_spectral_bounds::Bool=false,
    spectral_safety_margin::Float64=0.0,
    purification_method::Symbol=:sp2,
    chemical_potential::Union{Nothing,Real}=nothing,
    mcweeny_form::Symbol=:standard,
    mcweeny_trace_target::Union{Nothing,Real}=nothing,
    mcweeny_trace_tolerance::Union{Nothing,Real}=nothing,
    gc_policy::Symbol=:automatic,
    gc_period::Integer=10,
    gc_threshold_bytes::Integer=1 << 30,
    purification_cleanup::Function=() -> nothing,
    record_energy::Bool=false,
    stable_iterations::Integer=2,
    detect_two_cycles::Bool=true,
    two_cycle_tolerance::Union{Nothing,Real}=nothing,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    to_gpu=identity,
    to_cpu=identity,
    cleanup= () -> nothing)

    _validate_gc_policy(gc_policy, gc_period, gc_threshold_bytes)
    stable_iterations > 0 || throw(ArgumentError("stable_iterations must be positive"))


    if verbose == :all
        println(io, "Starting SCF iterations with parameters:")
        println(io, "  Max Iterations: $(sys.params.scf_max_iterations)")
        println(io, "  Convergence Tolerance: $(sys.params.scf_tol)")
        println(io, "  Mixing Parameter: $(sys.params.scf_mixing)")
    elseif verbose == :Density
        println(io, "Starting SCF iterations with parameters:")
        println(io, "  Max Iterations: $(sys.params.scf_max_iterations)")
        println(io, "  Convergence Tolerance: $(sys.params.scf_tol)")
        println(io, "  Mixing Parameter: $(sys.params.scf_mixing)")
    end
    verbose != :nothing && flush(io)

    sys.H0 = to_gpu(sys.H0)
    sys.VH = to_gpu(sys.VH)
    sys.VF = to_gpu(sys.VF) #

   
    params = sys.params
    residual_tolerance = params.scf_tol / 100
    cycle_tolerance = isnothing(two_cycle_tolerance) ? residual_tolerance : Float64(two_cycle_tolerance)
    isfinite(cycle_tolerance) && cycle_tolerance > 0 || throw(ArgumentError(
        "two_cycle_tolerance must be finite and positive",
    ))

    converged = false
    termination_reason = :max_iterations
    history = SCFIterationRecord[]
    sys.scf_diagnostics = SCFDiagnostics(history, false, :running)
    vh_residual = Inf
    vf_residual = Inf
    rho_residual = Inf
    commutator_residual = Inf
    rho_two_steps_ago = nothing
    overwrote_scf_progress = false
    mcweeny_identity_mpo = if purification_method == :mcweeny_mu &&
                              mcweeny_form == :horner
        to_gpu(Identity_MPO(sys.sites))
    else
        nothing
    end
    for iter in 1:params.scf_max_iterations
        # Step 1: Obtain density matrix!
        ρ0_device = construct_rho_0(
            sys, params, H_min, H_max;
            to_gpu=to_gpu,
            verify_spectral_bounds=verify_spectral_bounds,
            safety_margin=spectral_safety_margin,
            method=purification_method,
            chemical_potential=chemical_potential,
        )
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
            mcweeny_form=mcweeny_form,
            mcweeny_identity_mpo=mcweeny_identity_mpo,
            mcweeny_trace_target=mcweeny_trace_target,
            mcweeny_trace_tolerance=mcweeny_trace_tolerance,
            gc_policy=gc_policy,
            gc_period=gc_period,
            gc_threshold_bytes=gc_threshold_bytes,
            device_cleanup=purification_cleanup,
            io=io,
            overwrite_progress=overwrite_progress,
        )
        if !purification.converged && !allow_unconverged_purification
            push!(history, SCFIterationRecord(
                iter,
                purification.trace,
                vh_residual,
                vf_residual,
                rho_residual,
                commutator_residual,
                Inf,
                false,
                purification.termination_reason,
                purification.iterations,
                purification.selected_iteration,
                maxlinkdim(sys.ρ),
                maxlinkdim(sys.VH),
                maxlinkdim(sys.VF),
                nothing,
                nothing,
            ))
            sys.scf_diagnostics = SCFDiagnostics(history, false, :purification_failed)
            finish_iteration_progress(io, overwrote_scf_progress)
            verbose != :nothing && println(
                io, "SCF stopped: purification ended with $(purification.termination_reason).",
            )
            return false
        end
        ρ_purified_device = purification.rho

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
        vh_mpo_cpu, vf_mpo_cpu, fock_components = _extract_mean_fields_with_components(sys)

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
        two_cycle_residual = iter >= 3 ?
            _mpo_relative_change(sys.ρ, rho_two_steps_ago, params) : Inf
        energy_total = record_energy ? _scf_energy_total(sys) : nothing
        push!(history, SCFIterationRecord(
            iter,
            purification.trace,
            vh_residual,
            vf_residual,
            rho_residual,
            commutator_residual,
            two_cycle_residual,
            purification.converged,
            purification.termination_reason,
            purification.iterations,
            purification.selected_iteration,
            maxlinkdim(ρ_purified_device),
            maxlinkdim(vh_mpo),
            maxlinkdim(vf_mpo),
            maxlinkdim(H_effective),
            energy_total,
        ))
        if verbose != :nothing
            details = @sprintf(
                "Tr=%.12g | VH=%.3e | VF=%.3e | ρ=%.3e | [H,ρ]=%.3e",
                purification.trace,
                vh_residual,
                vf_residual,
                rho_residual,
                commutator_residual,
            )
            details *= @sprintf(
                " | χ=(ρ:%d,VH:%d,VF:%d,H:%d)",
                maxlinkdim(ρ_purified_device), maxlinkdim(vh_mpo),
                maxlinkdim(vf_mpo), maxlinkdim(H_effective),
            )
            if !isnothing(fock_components)
                horizontal_norm = sqrt(max(0.0, real(inner(fock_components.horizontal, fock_components.horizontal))))
                vertical_norm = sqrt(max(0.0, real(inner(fock_components.vertical, fock_components.vertical))))
                details *= @sprintf(" | ||VFₓ||=%.3e | ||VFᵧ||=%.3e", horizontal_norm, vertical_norm)
            end
            overwrote_scf_progress = print_iteration_progress(
                io, "SCF", iter, params.scf_max_iterations, details;
                overwrite=overwrite_progress,
            )
        end

        current_record = history[end]
        if detect_two_cycles && _scf_is_two_cycle(current_record, cycle_tolerance)
            termination_reason = :two_cycle_detected
            finish_iteration_progress(io, overwrote_scf_progress)
            verbose != :nothing && println(io, "SCF stopped: detected a two-cycle in the density matrix.")
            break
        end

        if _scf_has_stable_history(
            history,
            residual_tolerance;
            required_iterations=stable_iterations,
        )
            converged = true
            termination_reason = :converged
            finish_iteration_progress(io, overwrote_scf_progress)
            if verbose == :all
                println(io, "SCF converged in $iter iterations.")
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
        rho_two_steps_ago = rho_previous
        # Optional external cleanup may reclaim device memory; host GC follows
        # the explicit policy above instead of being forced every iteration.
        cleanup()
        maybe_collect_garbage!(iter;
            gc_policy=gc_policy,
            gc_period=gc_period,
            gc_threshold_bytes=gc_threshold_bytes,
        )
    end
    if !converged && verbose != :nothing
        if termination_reason == :max_iterations
            finish_iteration_progress(io, overwrote_scf_progress)
            println(io, "SCF did not converge within $(params.scf_max_iterations) iterations. Final maximum residual: $(maximum((vh_residual, vf_residual, rho_residual, commutator_residual)) * 100) %")
        end
    end
    sys.scf_diagnostics = SCFDiagnostics(history, converged, termination_reason)
    return converged    

end
