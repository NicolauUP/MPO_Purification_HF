function exact_occupied_projector(H, Ne)
    decomposition = eigen(Hermitian((H + H') / 2))
    occupied = decomposition.vectors[:, sortperm(decomposition.values)[1:Ne]]
    return occupied * occupied'
end

function exact_spectral_bounds(H; padding=0.5)
    values = eigvals(Hermitian((H + H') / 2))
    return (minimum(values) - padding, maximum(values) + padding)
end

function purify_static_system(sys; padding=0.5)
    H_eff = +(sys.H0, sys.VH, sys.VF;
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
    H_dense = dense_matrix(H_eff, sys)
    bounds = exact_spectral_bounds(H_dense; padding=padding)
    rho0 = construct_rho_0(
        sys, sys.params, bounds...; verify_spectral_bounds=true,
    )
    result = perform_purification(
        rho0, sys.params;
        verbose=0,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    return result, H_dense
end

@testset "M2.3 adaptive purification validation" begin
    diagonal_potential(x) = (-2.0, -0.5, 0.7, 2.0)[x + 1]
    for Ne in 1:3
        params = parameters_1d(
            t=0.0,
            W=diagonal_potential,
            U=0.0,
            density=Ne / 4,
            purification_steps=30,
        )
        sys = System(params)
        result, H_dense = purify_static_system(sys)
        rho_dense = dense_matrix(result.rho, sys)
        exact = exact_occupied_projector(H_dense, Ne)

        @test result.converged
        @test result.idempotency_residual < 1e-3
        @test isapprox(tr(rho_dense), Ne; atol=1e-6, rtol=1e-6)
        @test opnorm(rho_dense - rho_dense') < 1e-10
        @test minimum(eigvals(Hermitian(rho_dense))) >= -1e-10
        @test maximum(eigvals(Hermitian(rho_dense))) <= 1.0 + 1e-10
        @test opnorm(H_dense * rho_dense - rho_dense * H_dense) < 1e-10
        @test opnorm(rho_dense - exact) < 2e-3
    end

    shifted_scaled(x) = 1.7 * diagonal_potential(x) + 0.4
    scaled_params = parameters_1d(
        t=0.0,
        W=shifted_scaled,
        U=0.0,
        density=0.5,
        purification_steps=30,
    )
    scaled_sys = System(scaled_params)
    scaled_result, _ = purify_static_system(scaled_sys)
    reference_sys = System(parameters_1d(
        t=0.0,
        W=diagonal_potential,
        U=0.0,
        density=0.5,
        purification_steps=30,
    ))
    reference_result, _ = purify_static_system(reference_sys)
    @test isapprox(
        dense_matrix(scaled_result.rho, scaled_sys),
        dense_matrix(reference_result.rho, reference_sys);
        atol=2e-3,
        rtol=2e-3,
    )

    chain_params = parameters_1d(
        t=-0.7,
        W=x -> (0.2, -0.1, 0.05, 0.4)[x + 1],
        U=0.0,
        density=0.5,
        purification_steps=30,
    )
    chain_sys = System(chain_params)
    chain_result, chain_H = purify_static_system(chain_sys)
    chain_rho = dense_matrix(chain_result.rho, chain_sys)
    @test chain_result.converged
    @test chain_result.idempotency_residual < 1e-3
    @test opnorm(chain_H * chain_rho - chain_rho * chain_H) < 1e-8
    @test opnorm(chain_rho - exact_occupied_projector(chain_H, 2)) < 2e-3

    random_bonds = (-0.41, 0.18, -0.63, 0.29, -0.52, 0.07, 0.44)
    random_potential = (0.13, -0.31, 0.22, -0.08, 0.37, -0.19, 0.05, 0.28)
    random_params = parameters_1d(
        L=3,
        t=x -> x < 7 ? random_bonds[x + 1] : 0.0,
        W=x -> random_potential[x + 1],
        U=0.0,
        density=3 / 8,
        purification_steps=35,
        itensors_maxdim=64,
    )
    random_sys = System(random_params)
    random_result, random_H = purify_static_system(random_sys)
    random_rho = dense_matrix(random_result.rho, random_sys)
    @test random_result.converged
    @test random_result.idempotency_residual < 1e-3
    @test opnorm(random_H * random_rho - random_rho * random_H) < 1e-8
    @test opnorm(random_rho - exact_occupied_projector(random_H, 3)) < 2e-3

    square_params = parameters_square(
        t=(-0.6, -0.35),
        W=(x, y) -> 0.11x + 0.07y,
        U=0.0,
        density=0.5,
        purification_steps=35,
        itensors_maxdim=64,
    )
    square_sys = System(square_params)
    square_result, square_H = purify_static_system(square_sys)
    square_rho = dense_matrix(square_result.rho, square_sys)
    @test square_result.converged
    @test square_result.idempotency_residual < 1e-3
    @test opnorm(square_H * square_rho - square_rho * square_H) < 1e-7
    @test opnorm(square_rho - exact_occupied_projector(square_H, 8)) < 3e-3

    degenerate_params = parameters_1d(
        t=0.0,
        W=x -> (-1.0, 0.0, 0.0, 1.0)[x + 1],
        U=0.0,
        density=0.5,
        purification_steps=20,
    )
    degenerate_sys = System(degenerate_params)
    degenerate_H = dense_matrix(degenerate_sys.H0, degenerate_sys)
    degenerate_bounds = exact_spectral_bounds(degenerate_H)
    degenerate_rho0 = construct_rho_0(
        degenerate_sys, degenerate_params, degenerate_bounds...;
        verify_spectral_bounds=true,
    )
    degenerate = @test_logs (:warn, r"Purification did not converge") perform_purification(
        degenerate_rho0, degenerate_params;
        verbose=0,
        spectral_bounds=degenerate_bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    @test !degenerate.converged
    @test degenerate.termination_reason == :max_iterations

    near_degenerate_params = parameters_1d(
        t=0.0,
        W=x -> (-1.0, -1e-3, 1e-3, 1.0)[x + 1],
        U=0.0,
        density=0.5,
        purification_steps=10,
    )
    near_degenerate_sys = System(near_degenerate_params)
    near_degenerate_H = dense_matrix(near_degenerate_sys.H0, near_degenerate_sys)
    near_degenerate_bounds = exact_spectral_bounds(near_degenerate_H)
    near_degenerate_rho0 = construct_rho_0(
        near_degenerate_sys, near_degenerate_params, near_degenerate_bounds...;
        verify_spectral_bounds=true,
    )
    near_degenerate = @test_logs (:warn, r"Purification did not converge") perform_purification(
        near_degenerate_rho0, near_degenerate_params;
        verbose=0,
        spectral_bounds=near_degenerate_bounds,
        spectral_bounds_validation=:exact_small_system,
    )
    @test !near_degenerate.converged
    @test near_degenerate.termination_reason == :max_iterations
end
