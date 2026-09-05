function sp2_static_system()
    params = parameters_1d(
        t=0.0,
        W=x -> (-2.0, -0.5, 0.7, 2.0)[Int(x) + 1],
        U=0.0,
        density=0.5,
        purification_steps=30,
    )
    sys = System(params)
    H = dense_matrix(sys.H0, sys)
    bounds = exact_spectral_bounds(H; padding=0.5)
    rho0 = construct_rho_0(
        sys, params, bounds...;
        method=:sp2,
        verify_spectral_bounds=true,
    )
    return sys, params, H, bounds, rho0
end

@testset "M3.2 selectable SP2 backend" begin
    sys, params, H, bounds, rho0 = sp2_static_system()
    result = perform_purification(
        rho0, params;
        verbose=0,
        method=:sp2,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    rho = dense_matrix(result.rho, sys)
    exact = exact_occupied_projector(H, 2)

    @test result.method == :sp2
    @test result.converged
    @test result.termination_reason == :idempotency_threshold
    @test result.trace_error <= 1e-6
    @test result.idempotency_residual < 1e-3
    @test result.hermiticity_residual <= 1e-10
    @test result.work.squares == result.iterations
    @test result.work.cubes == 0
    @test result.work.max_bond_dimension >= result.final_bond_dimension
    @test length(result.work.bond_dimensions) == result.iterations
    @test all(dimensions -> dimensions[3] == 0, result.work.bond_dimensions)
    @test opnorm(rho - exact) < 2e-3

    strict_result = perform_purification(
        rho0, params;
        verbose=0,
        method=:sp2,
        spectral_bounds=bounds,
        sp2_idempotency_tolerance=1e-4,
    )
    @test strict_result.converged
    @test strict_result.idempotency_residual < 1e-4
    @test strict_result.iterations >= result.iterations

    relaxed_trace_result = perform_purification(
        rho0, params;
        verbose=0,
        method=:sp2,
        spectral_bounds=bounds,
        sp2_trace_tolerance=1e-2,
    )
    @test relaxed_trace_result.converged
    @test relaxed_trace_result.trace_error <= 1e-2
    @test relaxed_trace_result.iterations <= result.iterations

    default_result = perform_purification(
        rho0, params;
        verbose=0,
        spectral_bounds=bounds,
    )
    @test default_result.method == :sp2
    @test default_result.converged

    scf_sys = System(params)
    phase_records = NamedTuple[]
    detail_phase_records = NamedTuple[]
    synchronization_calls = Ref(0)
    @test run_scf!(
        scf_sys, bounds...;
        verbose=:nothing,
        verify_spectral_bounds=true,
        purification_method=:sp2,
        sp2_idempotency_tolerance=1e-4,
        sp2_trace_tolerance=1e-6,
        phase_callback=record -> push!(phase_records, record),
        detail_phase_callback=record -> push!(detail_phase_records, record),
        phase_synchronize=() -> (synchronization_calls[] += 1),
    )
    @test length(phase_records) == length(scf_diagnostics(scf_sys).history)
    @test synchronization_calls[] > 0
    @test all(record -> record.initialization_time_s >= 0, phase_records)
    @test all(record -> record.purification_time_s >= 0, phase_records)
    @test all(record -> record.density_to_host_time_s >= 0, phase_records)
    @test all(record -> record.mean_field_time_s >= 0, phase_records)
    @test all(record -> record.fields_to_device_time_s >= 0, phase_records)
    @test all(record -> record.device_diagnostics_time_s >= 0, phase_records)
    @test all(record -> record.residuals_time_s >= 0, phase_records)
    @test all(record -> record.energy_time_s >= 0, phase_records)
    @test all(record -> record.mixing_time_s >= 0, phase_records)
    @test all(record -> record.measured_iteration_time_s >= 0, phase_records)
    @test all(record -> record.purification_iterations > 0, phase_records)
    detail_phases = Set(record.phase for record in detail_phase_records)
    @test :square_apply in detail_phases
    @test :scalar_diagnostics in detail_phases
    @test :hartree in detail_phases
    @test :fock in detail_phases
    @test :effective_hamiltonian in detail_phases
    @test :commutator_Hrho in detail_phases
    @test :commutator_rhoH in detail_phases
    @test :commutator_difference in detail_phases
    @test :vh_residual in detail_phases
    @test :vf_residual in detail_phases
    @test :rho_residual in detail_phases
    @test all(record -> record.elapsed_time_s >= 0, detail_phase_records)
    @test opnorm(dense_matrix(scf_sys.ρ, scf_sys) - exact) < 2e-3
    returned_hamiltonian = dense_matrix(scf_sys.H0, scf_sys) +
        dense_matrix(scf_sys.VH, scf_sys) + dense_matrix(scf_sys.VF, scf_sys)
    returned_rho = dense_matrix(scf_sys.ρ, scf_sys)
    @test opnorm(returned_hamiltonian * returned_rho - returned_rho * returned_hamiltonian) < 1e-8

    @test_throws ArgumentError perform_purification(
        rho0, params; verbose=0, method=:sp2,
    )
    @test_throws ArgumentError perform_purification(
        rho0, params; verbose=0, method=:unsupported,
    )
    @test_throws ArgumentError perform_purification(
        rho0, params;
        verbose=0,
        method=:sp2,
        spectral_bounds=bounds,
        sp2_idempotency_tolerance=0.0,
    )
    @test_throws ArgumentError perform_purification(
        rho0, params;
        verbose=0,
        method=:sp2,
        spectral_bounds=bounds,
        sp2_trace_tolerance=-1.0,
    )
end
