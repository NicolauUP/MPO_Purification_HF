function dense_hf_energy_1d(H0, rho, U)
    N = size(H0, 1)
    interaction = U * sum(
        rho[i, i] * rho[i + 1, i + 1] - abs2(real(rho[i, i + 1]))
        for i in 1:(N - 1)
    )
    return real(tr(H0 * rho)) + interaction
end

@testset "M6.3 compact end-to-end 1D regressions" begin
    # Analytical N=2 reference: H0 = -(|1><2| + |2><1|), one particle.
    # The occupied bonding projector and total energy are both exact.
    for method in (:sp2, :palser_manolopoulos)
        params = parameters_1d(
            L=1,
            t=-1.0,
            U=0.0,
            density=0.5,
            purification_steps=25,
            scf_max_iterations=5,
        )
        sys = System(params)
        @test run_scf!(sys, -2.0, 2.0; purification_method=method)
        observables = observables_1d(sys)
        @test dense_matrix(sys.ρ, sys) ≈ [0.5 0.5; 0.5 0.5] atol=3e-3
        @test isapprox(observables.particle_number, 1.0; atol=3e-3)
        @test isapprox(observables.energy.total, -1.0; atol=3e-3)
        @test observables.hermiticity_residual < 1e-8
        @test observables.idempotency_residual < 3e-3
        @test observables.stationarity_residual < 1e-7
    end

    # Non-interacting modulated N=4 reference. Direct diagonalization is the
    # independent oracle for its projector and band energy.
    for method in (:sp2, :palser_manolopoulos)
        params = parameters_1d(
            t=x -> (-0.7, -0.3, -0.5)[Int(x) + 1],
            W=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
            U=0.0,
            purification_steps=30,
            scf_max_iterations=5,
        )
        sys = System(params)
        H0 = dense_matrix(sys.H0, sys)
        bounds = exact_spectral_bounds(H0; padding=0.5)
        exact = exact_occupied_projector(H0, 2)
        @test run_scf!(sys, bounds...; purification_method=method)
        observables = observables_1d(sys)
        @test opnorm(dense_matrix(sys.ρ, sys) - exact) < 3e-3
        @test isapprox(observables.energy.total, real(tr(H0 * exact)); atol=3e-3)
        @test observables.stationarity_residual < 1e-7
    end

    # Weakly interacting N=4 case: compare the complete SCF state and energy
    # with the independently implemented dense iteration in M5.4.
    for method in (:sp2, :palser_manolopoulos)
        params = parameters_1d(
            t=-0.7,
            W=x -> (-0.2, 0.1, -0.05, 0.25)[Int(x) + 1],
            U=0.3,
            S=nothing,
            purification_steps=35,
            scf_mixing=0.5,
            scf_tol=0.1,
            scf_max_iterations=20,
        )
        sys = System(params)
        H0 = dense_matrix(sys.H0, sys)
        rho_dense, _, _ = dense_hf_1d(H0, params.U, 2;
            mixing=params.scf_mixing,
            max_iterations=params.scf_max_iterations,
        )
        @test run_scf!(sys, -5.0, 5.0; purification_method=method)
        observables = observables_1d(sys)
        @test opnorm(dense_matrix(sys.ρ, sys) - rho_dense) < 4e-3
        @test isapprox(
            observables.energy.total,
            dense_hf_energy_1d(H0, rho_dense, params.U);
            atol=5e-3,
        )
        @test observables.stationarity_residual < 1e-3
    end
end
