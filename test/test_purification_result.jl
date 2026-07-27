function purification_initial_state(params)
    sys = System(params)
    rho0 = construct_rho_0(sys, params, -3.0, 3.0; method=:palser_manolopoulos)
    return sys, rho0
end

@testset "M2.1 purification result contract" begin
    successful_params = parameters_1d(U=0.0, purification_steps=20)
    _, successful_rho0 = purification_initial_state(successful_params)
    successful = perform_purification(successful_rho0, successful_params; method=:palser_manolopoulos, verbose=0)

    @test successful isa PurificationResult
    @test successful.method == :palser_manolopoulos
    legacy_alias = perform_purification(
        successful_rho0, successful_params;
        method=:adaptive_pm_mcweeny,
        verbose=0,
    )
    @test legacy_alias.method == :palser_manolopoulos
    @test successful.converged
    @test successful.termination_reason == :idempotency_threshold
    @test 1 <= successful.iterations <= successful_params.purification_steps
    @test successful.selected_iteration == successful.iterations
    @test successful.trace_error <= 1e-10
    @test successful.idempotency_residual <= 2e-3
    @test successful.hermiticity_residual <= 1e-12
    @test successful.final_bond_dimension == maxlinkdim(successful.rho)
    @test successful.spectral_bounds === nothing
    @test successful.spectral_bounds_validation == :not_provided
    @test successful.truncation_cutoff == successful_params.itensors_tol
    @test successful.maxdim == successful_params.itensors_maxdim
    @test successful.work.squares == successful.iterations
    @test successful.work.cubes == successful.iterations
    @test successful.work.max_bond_dimension >= successful.final_bond_dimension
    @test successful.work.mean_bond_dimension >= 1
    @test length(successful.work.bond_dimensions) == successful.iterations
    @test all(dimensions -> dimensions[3] > 0, successful.work.bond_dimensions)

    limited_params = parameters_1d(U=0.0, purification_steps=1)
    _, limited_rho0 = purification_initial_state(limited_params)
    limited = @test_logs (:warn, r"Purification did not converge") perform_purification(
        limited_rho0, limited_params; method=:palser_manolopoulos, verbose=0,
    )
    @test !limited.converged
    @test limited.termination_reason == :max_iterations
    @test limited.iterations == 1

    drift_params = parameters_1d(U=0.0)
    drift_sys = System(drift_params)
    drift = @test_logs (:warn, r"Trace has drifted") perform_purification(
        Identity_MPO(drift_sys.sites), drift_params; method=:palser_manolopoulos, verbose=0,
    )
    @test !drift.converged
    @test drift.termination_reason == :trace_drift
    @test isapprox(drift.trace_error, 2.0; atol=1e-12, rtol=1e-12)

    refusing_sys, _ = purification_initial_state(limited_params)
    @test_logs (:warn, r"Purification did not converge") begin
        @test !run_scf!(refusing_sys, -3.0, 3.0; purification_method=:palser_manolopoulos, verbose=:nothing)
    end
    @test norm(dense_matrix(refusing_sys.ρ, refusing_sys)) == 0.0

    diagnostic_params = parameters_1d(U=0.0, purification_steps=1, scf_max_iterations=1)
    diagnostic_sys, _ = purification_initial_state(diagnostic_params)
    @test_logs (:warn, r"Purification did not converge") begin
        @test !run_scf!(
            diagnostic_sys, -3.0, 3.0;
            verbose=:nothing,
            purification_method=:palser_manolopoulos,
            allow_unconverged_purification=true,
        )
    end
    @test norm(dense_matrix(diagnostic_sys.ρ, diagnostic_sys)) > 0.0
end
