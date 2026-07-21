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

function finite_gap_square_sp2_regression()
    params = parameters_square(
        t=(-0.6, -0.35),
        W=(x, y) -> 0.11x + 0.07y,
        U=0.0,
        density=0.5,
        purification_steps=35,
        itensors_maxdim=64,
    )
    H, pm, rho_pm, sp2, rho_sp2 = crossvalidate_purification(params)
    Ne = round(Int, size(H, 1) * params.density)
    return H, Ne, pm, rho_pm, sp2, rho_sp2
end

@testset "M3.4 PM/SP2 scientific cross-validation" begin
    cases = (
        ("1D diagonal, zero hopping", parameters_1d(
            t=0.0,
            W=x -> (-2.0, -0.5, 0.7, 2.0)[Int(x) + 1],
            U=0.0,
            density=0.5,
            purification_steps=30,
        )),
        ("2D square L=4, noninteracting (nondegenerate reference)", parameters_square(
            t=(-0.6, -0.35),
            # This deliberately matches the non-interacting square SCF reference.
            # The bilinear term lifts the accidental near-degeneracy at the
            # half-filled Fermi level, so an exact occupied projector is a
            # well-defined target for both canonical purification schemes.
            W=(x, y) -> 0.11x + 0.07y + 0.013x * y,
            U=0.0,
            density=0.5,
            purification_steps=35,
            itensors_maxdim=64,
        )),
        ("1D hopping, noninteracting", parameters_1d(
            t=-0.7,
            W=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
            U=0.0,
            density=0.5,
            purification_steps=30,
        )),
    )

    for (label, params) in cases
        println("M3.4 cross-validation case: ", label)
        flush(stdout)
        H, pm, rho_pm, sp2, rho_sp2 = crossvalidate_purification(params)
        println(
            "  PM:  converged=", pm.converged,
            " termination=", pm.termination_reason,
            " iterations=", pm.iterations,
            " trace_error=", pm.trace_error,
            " chi=", pm.final_bond_dimension,
        )
        println(
            "  SP2: converged=", sp2.converged,
            " termination=", sp2.termination_reason,
            " iterations=", sp2.iterations,
            " trace_error=", sp2.trace_error,
            " chi=", sp2.final_bond_dimension,
        )
        flush(stdout)
        Ne = round(Int, size(H, 1) * params.density)
        exact = exact_occupied_projector(H, Ne)

        @testset "$label" begin
            @test pm.converged && sp2.converged
            @test pm.idempotency_residual < 1e-3
            @test sp2.idempotency_residual < 1e-3
            @test isapprox(tr(rho_pm), Ne; atol=1e-6, rtol=1e-6)
            @test isapprox(tr(rho_sp2), Ne; atol=1e-6, rtol=1e-6)
            @test opnorm(rho_pm - exact) < 3e-3
            @test opnorm(rho_sp2 - exact) < 3e-3
            @test opnorm(rho_pm - rho_sp2) < 4e-3
            @test opnorm(H * rho_sp2 - rho_sp2 * H) < 1e-7
            @test isapprox(real(tr(H * rho_sp2)), real(tr(H * exact)); atol=3e-3, rtol=3e-3)
        end
    end

    @testset "2D square L=4, finite-gap SP2 MPO regression" begin
        H, Ne, pm, rho_pm, sp2, rho_sp2 = finite_gap_square_sp2_regression()
        spectrum = eigvals(Hermitian((H + H') / 2))
        fermi_gap = spectrum[Ne + 1] - spectrum[Ne]
        exact = exact_occupied_projector(H, Ne)

        # This is a finite-gap problem, rather than an intrinsically ambiguous
        # degenerate Fermi-level projector. PM remains the reference result.
        @test fermi_gap > 0.1
        @test pm.converged
        @test opnorm(rho_pm - exact) < 3e-3

        # Current MPO-SP2 stagnates despite the finite gap. Keep this visible as
        # an expected failure until the truncation/branch-stagnation mechanism
        # is corrected; an unexpected pass requires revisiting these checks.
        @test_broken sp2.converged
        @test_broken isapprox(tr(rho_sp2), Ne; atol=1e-6, rtol=1e-6)
        @test_broken opnorm(H * rho_sp2 - rho_sp2 * H) < 1e-7
    end

    @testset "2D square L=4, finite-gap SP2 tight-cutoff recovery" begin
        # The same affine finite-gap Hamiltonian as the expected-failure case
        # above. At maxdim=64 the default failure is therefore not a bond-cap
        # effect: reducing the MPO cutoff from 1e-12 to 1e-14 preserves the
        # final trace-correction steps and recovers the dense-SP2 projector.
        params = parameters_square(
            t=(-0.6, -0.35),
            W=(x, y) -> 0.11x + 0.07y,
            U=0.0,
            density=0.5,
            purification_steps=35,
            itensors_tol=1e-14,
            itensors_maxdim=64,
        )
        H, _, _, sp2, rho_sp2 = crossvalidate_purification(params)
        Ne = round(Int, size(H, 1) * params.density)
        exact = exact_occupied_projector(H, Ne)
        spectrum = eigvals(Hermitian((H + H') / 2))

        @test spectrum[Ne + 1] - spectrum[Ne] > 0.1
        @test sp2.converged
        @test sp2.termination_reason == :idempotency_threshold
        @test isapprox(tr(rho_sp2), Ne; atol=1e-6, rtol=1e-6)
        @test opnorm(rho_sp2 - exact) < 3e-3
        @test opnorm(H * rho_sp2 - rho_sp2 * H) < 1e-7
    end
end
