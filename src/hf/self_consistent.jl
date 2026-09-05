function _mpo_relative_change(candidate::MPO, reference::MPO, params::AbstractModelParameters)
    difference = +(candidate, -reference; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return safe_relative_change(
        sqrt(max(0.0, real(inner(difference, difference)))),
        sqrt(max(0.0, real(inner(reference, reference)))),
    )
end

function _scf_commutator_residual(
    H::MPO,
    rho::MPO,
    params::AbstractModelParameters;
    phase_callback::Union{Nothing,Function}=nothing,
    phase_synchronize::Function=() -> nothing,
)
    Hrho = _profiled_operation(
        :commutator_Hrho; callback=phase_callback, synchronize=phase_synchronize,
    ) do
        apply(H, rho; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end
    rhoH = _profiled_operation(
        :commutator_rhoH; callback=phase_callback, synchronize=phase_synchronize,
    ) do
        apply(rho, H; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end
    return _profiled_operation(
        :commutator_difference; callback=phase_callback, synchronize=phase_synchronize,
    ) do
        _mpo_relative_change(Hrho, rhoH, params)
    end
end

function _scf_energy_total(
    sys::System;
    hopping_hamiltonian::MPO=sys.H0,
)
    if sys.params isa Parameters1D
        return nearest_neighbor_hf_energy_1d(
            sys; hopping_hamiltonian=hopping_hamiltonian,
        ).total
    elseif sys.params isa ParametersSquare
        return nearest_neighbor_hf_energy_square(
            sys; hopping_hamiltonian=hopping_hamiltonian,
        ).total
    end
    return nothing
end

function _scf_record_within_tolerance(
    record::SCFIterationRecord,
    tolerance::Real;
    require_stationarity::Bool=true,
)
    residuals = (
        record.vh_residual,
        record.vf_residual,
        record.rho_residual,
    )
    require_stationarity && (residuals = (residuals..., record.commutator_residual))
    return all(isfinite, residuals) && maximum(residuals) <= tolerance
end

function _scf_has_stable_history(
    history::Vector{SCFIterationRecord},
    tolerance::Real;
    required_iterations::Integer=2,
    require_stationarity::Bool=true,
)
    required_iterations > 0 || throw(ArgumentError("required_iterations must be positive"))
    length(history) >= required_iterations || return false
    return all(record -> _scf_record_within_tolerance(
            record, tolerance; require_stationarity=require_stationarity,
        ),
        @view history[(end - required_iterations + 1):end])
end

function _scf_is_two_cycle(record::SCFIterationRecord, tolerance::Real)
    return isfinite(record.two_cycle_residual) &&
           record.two_cycle_residual <= tolerance &&
           isfinite(record.rho_residual) && record.rho_residual > tolerance
end

"""Safeguarded DIIS state for the coupled `(V_H, V_F)` MPO fixed-point map.

Each retained entry represents `F(V_H, V_F)` and its residual
`F(V_H, V_F) - (V_H, V_F)`. The vectors are kept on the host because the
field extractors already return host MPOs; only the accepted mixed fields are
transferred to the selected runtime for the following purification step.
"""
mutable struct MPOPulayMixer
    outputs::Vector{Tuple{MPO,MPO}}
    residuals::Vector{Tuple{MPO,MPO}}
    history::Int
    warmup::Int
    regularization::Float64
    coefficient_limit::Float64
    step_limit::Float64
end

function MPOPulayMixer(; history::Integer=4, warmup::Integer=4,
    regularization::Real=1e-12, coefficient_limit::Real=8.0,
    step_limit::Real=20.0)
    history >= 2 || throw(ArgumentError("Pulay history must be at least 2"))
    warmup >= 2 || throw(ArgumentError("Pulay warmup must be at least 2"))
    warmup <= history || throw(ArgumentError("Pulay warmup must not exceed history"))
    isfinite(regularization) && regularization >= 0 || throw(ArgumentError(
        "Pulay regularization must be finite and nonnegative",
    ))
    isfinite(coefficient_limit) && coefficient_limit > 0 || throw(ArgumentError(
        "Pulay coefficient limit must be finite and positive",
    ))
    isfinite(step_limit) && step_limit > 0 || throw(ArgumentError(
        "Pulay step limit must be finite and positive",
    ))
    return MPOPulayMixer(
        Tuple{MPO,MPO}[], Tuple{MPO,MPO}[], Int(history), Int(warmup),
        Float64(regularization), Float64(coefficient_limit), Float64(step_limit),
    )
end

function _mpo_field_difference(candidate::Tuple{MPO,MPO}, reference::Tuple{MPO,MPO}, params)
    return (
        +(candidate[1], -reference[1]; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim),
        +(candidate[2], -reference[2]; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim),
    )
end

_mpo_field_inner(left::Tuple{MPO,MPO}, right::Tuple{MPO,MPO}) =
    real(inner(left[1], right[1]) + inner(left[2], right[2]))

_mpo_field_norm(fields::Tuple{MPO,MPO}) = sqrt(max(0.0, _mpo_field_inner(fields, fields)))

function _mpo_field_linear_combination(
    input::Tuple{MPO,MPO}, outputs::Vector{Tuple{MPO,MPO}}, weights::AbstractVector,
    damping::Real, params,
)
    hartree = (1 - damping) * input[1]
    fock = (1 - damping) * input[2]
    for index in eachindex(outputs, weights)
        hartree = +(hartree, damping * weights[index] * outputs[index][1];
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        fock = +(fock, damping * weights[index] * outputs[index][2];
            cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end
    return (hartree, fock)
end

"""Return a safe DIIS field update, or a linear fallback and its reason."""
function _pulay_update!(
    mixer::MPOPulayMixer, input::Tuple{MPO,MPO}, output::Tuple{MPO,MPO}, params;
    damping::Real,
)
    isfinite(damping) && 0 < damping <= 1 || throw(ArgumentError(
        "Pulay damping must be finite and lie in (0, 1]",
    ))
    residual = _mpo_field_difference(output, input, params)
    push!(mixer.outputs, output)
    push!(mixer.residuals, residual)
    while length(mixer.outputs) > mixer.history
        popfirst!(mixer.outputs)
        popfirst!(mixer.residuals)
    end

    count = length(mixer.residuals)
    linear = _mpo_field_linear_combination(input, [output], [1.0], damping, params)
    count < mixer.warmup && return linear, :linear_warmup

    gram = Matrix{Float64}(undef, count, count)
    for row in 1:count, column in 1:count
        gram[row, column] = _mpo_field_inner(mixer.residuals[row], mixer.residuals[column])
    end
    all(isfinite, gram) || return linear, :linear_nonfinite_gram_fallback
    scale = max(maximum(abs, gram), 1.0)
    @inbounds for index in 1:count
        gram[index, index] += mixer.regularization * scale
    end
    system = zeros(Float64, count + 1, count + 1)
    system[1:count, 1:count] .= gram
    system[1:count, count + 1] .= 1.0
    system[count + 1, 1:count] .= 1.0
    coefficients = try
        system \ vcat(zeros(Float64, count), 1.0)
    catch
        return linear, :linear_solve_fallback
    end
    weights = coefficients[1:count]
    all(isfinite, weights) || return linear, :linear_nonfinite_coefficient_fallback
    maximum(abs, weights) <= mixer.coefficient_limit ||
        return linear, :linear_coefficient_fallback

    candidate = _mpo_field_linear_combination(
        input, mixer.outputs, weights, damping, params,
    )
    candidate_step = _mpo_field_norm(_mpo_field_difference(candidate, input, params))
    residual_norm = _mpo_field_norm(residual)
    isfinite(candidate_step) && candidate_step <= mixer.step_limit * max(residual_norm, eps(Float64)) ||
        return linear, :linear_step_fallback
    return candidate, :pulay
end

function run_scf!(sys::System, H_min::Float64, H_max::Float64; 
    verbose::Symbol=:nothing,
    allow_unconverged_purification::Bool=false,
    verify_spectral_bounds::Bool=false,
    spectral_safety_margin::Float64=0.0,
    purification_method::Symbol=:sp2,
    square_fock_method::Symbol=:binary_carry,
    sp2_idempotency_tolerance::Real=1e-3,
    sp2_trace_tolerance::Union{Nothing,Real}=nothing,
    chemical_potential::Union{Nothing,Real}=nothing,
    mcweeny_form::Symbol=:standard,
    mcweeny_trace_target::Union{Nothing,Real}=nothing,
    mcweeny_trace_tolerance::Union{Nothing,Real}=nothing,
    gc_policy::Symbol=:automatic,
    gc_period::Integer=10,
    gc_threshold_bytes::Integer=1 << 30,
    purification_cleanup::Function=() -> nothing,
    record_energy::Bool=false,
    mixing_method::Symbol=:linear,
    pulay_history::Integer=4,
    pulay_warmup::Integer=4,
    pulay_regularization::Real=1e-12,
    pulay_coefficient_limit::Real=8.0,
    pulay_step_limit::Real=20.0,
    stable_iterations::Integer=2,
    require_stationarity::Bool=true,
    measure_stationarity::Bool=true,
    detect_two_cycles::Bool=true,
    two_cycle_tolerance::Union{Nothing,Real}=nothing,
    io::IO=stdout,
    overwrite_progress::Bool=io isa Base.TTY,
    to_gpu=identity,
    to_cpu=identity,
    phase_callback::Union{Nothing,Function}=nothing,
    detail_phase_callback::Union{Nothing,Function}=nothing,
    phase_synchronize::Function=() -> nothing,
    sp2_snapshot_callback::Union{Nothing,Function}=nothing,
    cleanup= () -> nothing)

    _validate_gc_policy(gc_policy, gc_period, gc_threshold_bytes)
    stable_iterations > 0 || throw(ArgumentError("stable_iterations must be positive"))
    require_stationarity && !measure_stationarity && throw(ArgumentError(
        "measure_stationarity=false is incompatible with require_stationarity=true",
    ))
    square_fock_method in (:tci, :binary_carry) || throw(ArgumentError(
        "square_fock_method must be :tci or :binary_carry, got $square_fock_method",
    ))
    mixing_method in (:linear, :pulay) || throw(ArgumentError(
        "mixing_method must be :linear or :pulay, got $mixing_method",
    ))


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

    # The extraction kernels produce fields on the host. Retaining the host
    # input pair allows Pulay to form its residual Gram matrix without extra
    # device-to-host transfers; the accepted fields are sent to the runtime
    # below for purification.
    host_fields = (sys.VH, sys.VF)
    pulay_mixer = mixing_method == :pulay ? MPOPulayMixer(
        history=pulay_history, warmup=pulay_warmup,
        regularization=pulay_regularization,
        coefficient_limit=pulay_coefficient_limit,
        step_limit=pulay_step_limit,
    ) : nothing

    # The density is explicitly returned to the host before Hartree/Fock
    # extraction. Keep the matching host H0 for optional per-iteration energy
    # accounting; otherwise record_energy would contract a CUDA H0 with a host
    # density MPO.
    H0_host = sys.H0
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
    timed_phase = function (operation)
        if isnothing(phase_callback)
            return operation(), NaN
        end
        phase_synchronize()
        started = time_ns()
        value = operation()
        phase_synchronize()
        return value, (time_ns() - started) / 1e9
    end
    for iter in 1:params.scf_max_iterations
        iteration_started = time_ns()
        # Step 1: Obtain density matrix!
        ρ0_device, initialization_time = timed_phase() do
            construct_rho_0(
                sys, params, H_min, H_max;
                to_gpu=to_gpu,
                verify_spectral_bounds=verify_spectral_bounds,
                safety_margin=spectral_safety_margin,
                method=purification_method,
                chemical_potential=chemical_potential,
            )
        end
        purification_verbose = verbose == :nothing ? 0 : 1
        purification, purification_time = timed_phase() do
            perform_purification(
                ρ0_device, params;
                verbose=purification_verbose,
                spectral_bounds=(H_min, H_max),
                spectral_bounds_validation=(
                    verify_spectral_bounds ? :exact_small_system : :supplied_unverified
                ),
                method=purification_method,
                sp2_idempotency_tolerance=sp2_idempotency_tolerance,
                sp2_trace_tolerance=sp2_trace_tolerance,
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
                phase_callback=isnothing(detail_phase_callback) ? nothing : event ->
                    detail_phase_callback(merge((scf_iteration=iter,), event)),
                phase_synchronize=phase_synchronize,
                sp2_snapshot_callback=isnothing(sp2_snapshot_callback) ? nothing : event ->
                    sp2_snapshot_callback(merge((scf_iteration=iter,), event)),
            )
        end
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
        rho_host, density_to_host_time = timed_phase(() -> to_cpu(ρ_purified_device))
        sys.ρ = rho_host

        # Step 2: Extract Hartree potential
        mean_fields, mean_field_time = timed_phase(
            () -> _extract_mean_fields_with_components(
                sys; square_fock_method=square_fock_method,
                phase_callback=isnothing(detail_phase_callback) ? nothing : event ->
                    detail_phase_callback(merge((scf_iteration=iter,), event)),
                phase_synchronize=phase_synchronize,
            ),
        )
        vh_mpo_cpu, vf_mpo_cpu, fock_components = mean_fields

        device_fields, fields_to_device_time = timed_phase() do
            (to_gpu(vh_mpo_cpu), to_gpu(vf_mpo_cpu))
        end
        vh_mpo, vf_mpo = device_fields


        # Step 4: Check convergence
        device_diagnostics, device_diagnostics_time = timed_phase() do
            H = _profiled_operation(
                :effective_hamiltonian;
                callback=detail_phase_callback,
                synchronize=phase_synchronize,
                metadata=(scf_iteration=iter,),
            ) do
                +(sys.H0, sys.VH, sys.VF;
                    cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
                )
            end
            residual = measure_stationarity ?
                _scf_commutator_residual(
                    H, ρ_purified_device, params;
                    phase_callback=isnothing(detail_phase_callback) ? nothing : event ->
                        detail_phase_callback(merge((scf_iteration=iter,), event)),
                    phase_synchronize=phase_synchronize,
                ) : NaN
            (H, residual)
        end
        H_effective, commutator_residual = device_diagnostics
        residuals, residuals_time = timed_phase() do
            vh_change = _profiled_operation(
                :vh_residual; callback=detail_phase_callback,
                synchronize=phase_synchronize, metadata=(scf_iteration=iter,),
            ) do
                iter > 1 ? _mpo_relative_change(vh_mpo, sys.VH, params) : Inf
            end
            vf_change = _profiled_operation(
                :vf_residual; callback=detail_phase_callback,
                synchronize=phase_synchronize, metadata=(scf_iteration=iter,),
            ) do
                iter > 1 ? _mpo_relative_change(vf_mpo, sys.VF, params) : Inf
            end
            rho_change = _profiled_operation(
                :rho_residual; callback=detail_phase_callback,
                synchronize=phase_synchronize, metadata=(scf_iteration=iter,),
            ) do
                iter > 1 ? _mpo_relative_change(sys.ρ, rho_previous, params) : Inf
            end
            two_cycle_change = _profiled_operation(
                :two_cycle_residual; callback=detail_phase_callback,
                synchronize=phase_synchronize, metadata=(scf_iteration=iter,),
            ) do
                iter >= 3 ? _mpo_relative_change(sys.ρ, rho_two_steps_ago, params) : Inf
            end
            (vh_change, vf_change, rho_change, two_cycle_change)
        end
        vh_residual, vf_residual, rho_residual, two_cycle_residual = residuals
        energy_total, energy_time = timed_phase() do
            if record_energy
                _profiled_operation(
                    :energy; callback=detail_phase_callback,
                    synchronize=phase_synchronize, metadata=(scf_iteration=iter,),
                ) do
                    _scf_energy_total(sys; hopping_hamiltonian=H0_host)
                end
            else
                nothing
            end
        end
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
        two_cycle_detected = detect_two_cycles &&
            _scf_is_two_cycle(current_record, cycle_tolerance)
        stable_history = !two_cycle_detected && _scf_has_stable_history(
            history,
            residual_tolerance;
            required_iterations=stable_iterations,
            require_stationarity=require_stationarity,
        )
        mixing_time = 0.0
        applied_mixing_method = :not_applied
        if !two_cycle_detected && !stable_history
            _, mixing_time = timed_phase() do
                if iter == 1
                    initial_host_fields = host_fields
                    sys.VH = vh_mpo
                    sys.VF = vf_mpo
                    host_fields = (vh_mpo_cpu, vf_mpo_cpu)
                    # Seed the DIIS history with the direct first update,
                    # preserving the established first-iteration behavior.
                    if !isnothing(pulay_mixer)
                        _pulay_update!(
                            pulay_mixer, initial_host_fields, host_fields, params;
                            damping=params.scf_mixing,
                        )
                    end
                    applied_mixing_method = :direct
                elseif isnothing(pulay_mixer)
                    sys.VH = +(
                        sys.params.scf_mixing * vh_mpo,
                        (1 - params.scf_mixing) * sys.VH;
                        cutoff=sys.params.itensors_tol,
                        maxdim=sys.params.itensors_maxdim,
                    )
                    sys.VF = +(
                        sys.params.scf_mixing * vf_mpo,
                        (1 - params.scf_mixing) * sys.VF;
                        cutoff=sys.params.itensors_tol,
                        maxdim=sys.params.itensors_maxdim,
                    )
                    applied_mixing_method = :linear
                else
                    mixed_fields, method = _pulay_update!(
                        pulay_mixer, host_fields, (vh_mpo_cpu, vf_mpo_cpu), params;
                        damping=params.scf_mixing,
                    )
                    sys.VH, sys.VF = to_gpu(mixed_fields[1]), to_gpu(mixed_fields[2])
                    host_fields = mixed_fields
                    applied_mixing_method = method
                end
            end
            rho_two_steps_ago = rho_previous
        end
        if !isnothing(phase_callback)
            phase_synchronize()
            phase_callback((
                iteration=iter,
                initialization_time_s=initialization_time,
                purification_time_s=purification_time,
                density_to_host_time_s=density_to_host_time,
                mean_field_time_s=mean_field_time,
                fields_to_device_time_s=fields_to_device_time,
                device_diagnostics_time_s=device_diagnostics_time,
                residuals_time_s=residuals_time,
                energy_time_s=energy_time,
                mixing_time_s=mixing_time,
                measured_iteration_time_s=(time_ns() - iteration_started) / 1e9,
                purification_iterations=purification.iterations,
                rho_bond_dimension=maxlinkdim(ρ_purified_device),
                hartree_bond_dimension=maxlinkdim(vh_mpo),
                fock_bond_dimension=maxlinkdim(vf_mpo),
                effective_hamiltonian_bond_dimension=maxlinkdim(H_effective),
                mixing_method=applied_mixing_method,
            ))
        end

        if two_cycle_detected
            termination_reason = :two_cycle_detected
            finish_iteration_progress(io, overwrote_scf_progress)
            verbose != :nothing && println(io, "SCF stopped: detected a two-cycle in the density matrix.")
            break
        end

        if stable_history
            converged = true
            termination_reason = :converged
            finish_iteration_progress(io, overwrote_scf_progress)
            if verbose == :all
                println(io, "SCF converged in $iter iterations.")
            end
            break
        end

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
