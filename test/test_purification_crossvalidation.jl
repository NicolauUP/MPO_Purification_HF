function crossvalidate_purification(params)
    pm_sys = System(params)
    H = dense_matrix(pm_sys.H0, pm_sys)
    bounds = exact_spectral_bounds(H; padding=0.5)
    pm_rho0 = construct_rho_0(
        pm_sys, params, bounds...;
        method=:palser_manolopoulos,
        verify_spectral_bounds=true,
    )
    pm = perform_purification(
        pm_rho0, params;
        method=:palser_manolopoulos,
        verbose=0,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )

    sp2_sys = System(params)
    sp2_rho0 = construct_rho_0(
        sp2_sys, params, bounds...;
        method=:sp2,
        verify_spectral_bounds=true,
    )
    sp2 = perform_purification(
        sp2_rho0, params;
        method=:sp2,
        verbose=0,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    return H, pm, dense_matrix(pm.rho, pm_sys), sp2, dense_matrix(sp2.rho, sp2_sys)
end

@testset "M3.4 PM/SP2 scientific cross-validation" begin
    cases = (
        parameters_1d(
            t=0.0,
            W=x -> (-2.0, -0.5, 0.7, 2.0)[Int(x) + 1],
            U=0.0,
            density=0.5,
            purification_steps=30,
        ),
        parameters_1d(
            t=-0.7,
            W=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
            U=0.0,
            density=0.5,
            purification_steps=30,
        ),
    )

    for params in cases
        H, pm, rho_pm, sp2, rho_sp2 = crossvalidate_purification(params)
        exact = exact_occupied_projector(H, 2)

        @test pm.converged && sp2.converged
        @test pm.idempotency_residual < 1e-3
        @test sp2.idempotency_residual < 1e-3
        @test isapprox(tr(rho_pm), 2.0; atol=1e-6, rtol=1e-6)
        @test isapprox(tr(rho_sp2), 2.0; atol=1e-6, rtol=1e-6)
        @test opnorm(rho_pm - exact) < 2e-3
        @test opnorm(rho_sp2 - exact) < 2e-3
        @test opnorm(rho_pm - rho_sp2) < 3e-3
        @test opnorm(H * rho_sp2 - rho_sp2 * H) < 1e-8
        @test isapprox(real(tr(H * rho_sp2)), real(tr(H * exact)); atol=2e-3, rtol=2e-3)
    end
end
