@testset "McWeeny standard and Horner polynomial forms" begin
    dense_rho = [
        0.25 0.03 0.00 0.01
        0.03 0.45 0.02 0.00
        0.00 0.02 0.65 0.04
        0.01 0.00 0.04 0.85
    ]
    dense_square = dense_rho * dense_rho
    dense_standard = 3dense_square - 2dense_square * dense_rho
    dense_horner = dense_square * (3I - 2dense_rho)
    @test dense_standard ≈ dense_horner atol=1e-14 rtol=1e-14

    params = parameters_1d(
        t=0.0,
        W=x -> (-1.5, -0.4, 0.6, 1.4)[Int(x) + 1],
        U=0.0,
        density=0.5,
        itensors_tol=1e-14,
        itensors_maxdim=64,
    )
    sys = System(params)
    rho = construct_rho_0(
        sys, params, -2.0, 2.0;
        method=:mcweeny_mu,
        chemical_potential=0.0,
    )
    rho_squared = apply(rho, rho;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    standard = MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:standard,
    )
    horner = MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params;
        form=:horner,
        identity_mpo=Identity_MPO(sys.sites),
    )

    @test standard.auxiliary_kind == :rho_cubed
    @test horner.auxiliary_kind == :horner_factor
    @test opnorm(dense_matrix(standard.rho, sys) - dense_matrix(horner.rho, sys)) < 1e-11
    @test_throws ArgumentError MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:horner,
    )
    @test_throws ArgumentError MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:unknown,
    )

    horner_result = perform_purification(
        rho, params;
        method=:mcweeny_mu,
        chemical_potential=0.0,
        spectral_bounds=(-2.0, 2.0),
        mcweeny_form=:horner,
        mcweeny_identity_mpo=Identity_MPO(sys.sites),
        verbose=0,
    )
    @test horner_result.converged
    @test horner_result.idempotency_residual < 1e-6

    scf_sys = System(params)
    @test run_scf!(
        scf_sys, -2.0, 2.0;
        purification_method=:mcweeny_mu,
        chemical_potential=0.0,
        mcweeny_form=:horner,
        verbose=:nothing,
    )
end

@testset "McWeeny retains the best trace-eligible capped iterate" begin
    params = parameters_square(
        L=4,
        t=(-0.6, -0.35),
        U=0.0,
        W=nothing,
        S=(x, y) -> iseven(Int(x) + Int(y)) ? 0.1 : -0.1,
        density=0.5,
        purification_steps=20,
        itensors_tol=1e-8,
        itensors_maxdim=16,
    )
    sys = System(params)
    rho0 = construct_rho_0(
        sys, params, -3.0, 3.0;
        method=:mcweeny_mu,
        chemical_potential=0.0,
    )

    unconstrained = @test_logs (:warn, r"returning selected iteration") perform_purification(
        rho0, params;
        method=:mcweeny_mu,
        chemical_potential=0.0,
        spectral_bounds=(-3.0, 3.0),
        verbose=0,
    )
    @test !unconstrained.converged
    @test unconstrained.termination_reason == :best_iterate_after_max_iterations
    @test unconstrained.iterations == 20
    @test unconstrained.selected_iteration == 13
    @test unconstrained.idempotency_residual < 2e-4
    @test unconstrained.work.squares == 20
    @test unconstrained.work.cubes == 19

    guarded = @test_logs (:warn, r"returning selected iteration") perform_purification(
        rho0, params;
        method=:mcweeny_mu,
        chemical_potential=0.0,
        spectral_bounds=(-3.0, 3.0),
        mcweeny_trace_target=8.0,
        mcweeny_trace_tolerance=1e-5,
        verbose=0,
    )
    @test !guarded.converged
    @test guarded.selected_iteration == 10
    @test guarded.target_particles == 8.0
    @test guarded.trace_error <= 1e-5
    @test guarded.idempotency_residual < 1e-3

    @test_throws ArgumentError perform_purification(
        rho0, params;
        method=:mcweeny_mu,
        chemical_potential=0.0,
        spectral_bounds=(-3.0, 3.0),
        mcweeny_trace_target=8.0,
        verbose=0,
    )
    @test_throws ArgumentError perform_purification(
        rho0, params;
        method=:mcweeny_mu,
        chemical_potential=0.0,
        spectral_bounds=(-3.0, 3.0),
        mcweeny_trace_tolerance=1e-5,
        verbose=0,
    )
end
