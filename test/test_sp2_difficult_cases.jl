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

    _, _, degenerate = sp2_diagonal_result((-1.0, 0.0, 0.0, 1.0), 2; steps=20)
    @test !degenerate.converged
    @test degenerate.termination_reason in (:stagnation, :max_iterations)

    _, _, limited = sp2_diagonal_result(values, 2; steps=1)
    @test !limited.converged
    @test limited.termination_reason == :max_iterations
end
