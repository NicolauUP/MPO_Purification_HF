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

    default_result = perform_purification(
        rho0, params;
        verbose=0,
        spectral_bounds=bounds,
    )
    @test default_result.method == :sp2
    @test default_result.converged

    scf_sys = System(params)
    @test run_scf!(
        scf_sys, bounds...;
        verbose=:nothing,
        verify_spectral_bounds=true,
        purification_method=:sp2,
    )
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
end
