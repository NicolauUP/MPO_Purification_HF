function sp2_diagonal_result(values, Ne; steps=35)
    params = parameters_1d(
        t=0.0,
        W=x -> values[Int(x) + 1],
        U=0.0,
        density=Ne / length(values),
        purification_steps=steps,
    )
    sys = System(params)
    H = dense_matrix(sys.H0, sys)
    bounds = exact_spectral_bounds(H; padding=0.5)
    rho0 = construct_rho_0(sys, params, bounds...; method=:sp2, verify_spectral_bounds=true)
    result = perform_purification(
        rho0, params;
        method=:sp2,
        verbose=0,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    return sys, H, result
end

@testset "M3.5 difficult SP2 cases" begin
    values = (-2.0, -0.5, 0.7, 2.0)
    for Ne in (1, 3)
        sys, H, result = sp2_diagonal_result(values, Ne)
        rho = dense_matrix(result.rho, sys)
        @test result.converged
        @test isapprox(tr(rho), Ne; atol=1e-6, rtol=1e-6)
        @test opnorm(rho - exact_occupied_projector(H, Ne)) < 2e-3
    end

    degenerate_result = @test_logs (:warn, r"SP2 purification did not converge") begin
        sp2_diagonal_result((-1.0, 0.0, 0.0, 1.0), 2; steps=20)
    end
    _, _, degenerate = degenerate_result
    @test !degenerate.converged
    @test degenerate.termination_reason in (:stagnation, :max_iterations)

    limited_result = @test_logs (:warn, r"SP2 purification did not converge") begin
        sp2_diagonal_result(values, 2; steps=1)
    end
    limited_sys, limited_H, limited = limited_result
    @test !limited.converged
    @test limited.termination_reason == :max_iterations
    @test limited.selected_iteration == 1
    limited_bounds = exact_spectral_bounds(limited_H; padding=0.5)
    limited_rho0 = construct_rho_0(
        limited_sys, limited_sys.params, limited_bounds...;
        method=:sp2,
        verify_spectral_bounds=true,
    )
    @test opnorm(
        dense_matrix(limited.rho, limited_sys) -
        dense_matrix(limited_rho0, limited_sys),
    ) < 1e-12

    # Same width, distinct positive Fermi gaps: both must select the exact rank-2 projector.
    # The 2e-3 gap needs 45 SP2 polynomials; it is finite, unlike the exactly
    # degenerate case above, so it must converge when given that budget.
    for (gap_values, steps) in (
        ((-2.0, -1.0, 1.0, 2.0), 35),
        ((-2.0, -1e-3, 1e-3, 2.0), 45),
    )
        sys, H, result = sp2_diagonal_result(gap_values, 2; steps)
        @test result.converged
        @test opnorm(dense_matrix(result.rho, sys) - exact_occupied_projector(H, 2)) < 2e-3
    end

    # H -> aH + bI with transformed bounds must preserve the occupied projector.
    reference_sys, _, reference = sp2_diagonal_result(values, 2)
    shifted_scaled = tuple((1.7 * value + 0.4 for value in values)...)
    transformed_sys, _, transformed = sp2_diagonal_result(shifted_scaled, 2)
    @test reference.converged && transformed.converged
    @test opnorm(
        dense_matrix(reference.rho, reference_sys) - dense_matrix(transformed.rho, transformed_sys),
    ) < 2e-3
end
