@testset "M2.2 spectral-bound contract" begin
    params = parameters_1d(U=0.0, purification_steps=20)
    sys = System(params)
    H_eff = +(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    dense_H = dense_matrix(H_eff, sys)
    eigenvalues = eigvals(Hermitian(dense_H))
    lambda_min, lambda_max = extrema(eigenvalues)
    bounds = (lambda_min - 0.25, lambda_max + 0.25)

    @test validate_spectral_bounds(bounds...) == bounds
    @test_throws ArgumentError validate_spectral_bounds(1.0, 1.0)
    @test_throws ArgumentError validate_spectral_bounds(NaN, 1.0)
    @test_throws ArgumentError validate_spectral_bounds(-1.0, Inf)
    @test_throws ArgumentError validate_spectral_bounds(-1.0, 1.0; safety_margin=-1e-3)

    exact_extrema = verify_spectral_bounds_exact(sys, H_eff, bounds...)
    @test isapprox(exact_extrema[1], lambda_min; atol=1e-12, rtol=1e-12)
    @test isapprox(exact_extrema[2], lambda_max; atol=1e-12, rtol=1e-12)
    @test_throws ArgumentError verify_spectral_bounds_exact(
        sys, H_eff, lambda_min + 0.1, lambda_max + 0.25,
    )
    @test_throws ArgumentError verify_spectral_bounds_exact(
        sys, H_eff, lambda_min - 0.25, lambda_max - 0.1,
    )
    @test_throws ArgumentError verify_spectral_bounds_exact(
        sys, H_eff, bounds...; safety_margin=0.5,
    )

    rho0 = construct_rho_0(
        sys, params, bounds...; verify_spectral_bounds=true,
    )
    scaled_eigenvalues = eigvals(Hermitian(dense_matrix(rho0, sys)))
    @test minimum(scaled_eigenvalues) >= -1e-12
    @test maximum(scaled_eigenvalues) <= 1.0 + 1e-12
    @test_throws ArgumentError construct_rho_0(
        sys, params, lambda_min + 0.1, lambda_max + 0.25;
        verify_spectral_bounds=true,
    )

    result = perform_purification(
        rho0, params;
        verbose=0,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    @test result.spectral_bounds == bounds
    @test result.spectral_bounds_validation == :exact_small_system
    @test_throws ArgumentError perform_purification(
        rho0, params; verbose=0, spectral_bounds=(1.0, 1.0),
    )

    scf_sys = System(parameters_1d(U=0.0, scf_max_iterations=3))
    @test run_scf!(
        scf_sys, bounds...;
        verbose=:nothing,
        verify_spectral_bounds=true,
    )
end
