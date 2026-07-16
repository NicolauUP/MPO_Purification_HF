function mcweeny_mu_system(; shift=0.0)
    params = parameters_1d(
        t=0.0,
        W=x -> (-2.0, -0.5, 0.7, 2.0)[Int(x) + 1] + shift,
        U=0.0,
        density=0.5,
        purification_steps=30,
    )
    sys = System(params)
    return sys, params, dense_matrix(sys.H0, sys)
end

function mcweeny_mu_result(sys, params, bounds, mu)
    rho0 = construct_rho_0(
        sys, params, bounds...;
        method=:mcweeny_mu,
        chemical_potential=mu,
        verify_spectral_bounds=true,
    )
    return perform_purification(
        rho0, params;
        method=:mcweeny_mu,
        chemical_potential=mu,
        spectral_bounds=bounds,
        spectral_bounds_validation=:exact_small_system,
        verbose=0,
    )
end

@testset "M3.2a direct McWeeny chemical-potential purification" begin
    sys, params, H = mcweeny_mu_system()
    bounds = exact_spectral_bounds(H; padding=0.5)

    half_filled = mcweeny_mu_result(sys, params, bounds, 0.0)
    rho_half = dense_matrix(half_filled.rho, sys)
    @test half_filled.converged
    @test half_filled.method == :mcweeny_mu
    @test half_filled.target_particles === nothing
    @test half_filled.trace_error === nothing
    @test isapprox(tr(rho_half), 2.0; atol=1e-6, rtol=1e-6)
    @test opnorm(rho_half - exact_occupied_projector(H, 2)) < 2e-3

    one_particle = mcweeny_mu_result(sys, params, bounds, -1.0)
    rho_one = dense_matrix(one_particle.rho, sys)
    @test one_particle.converged
    @test isapprox(tr(rho_one), 1.0; atol=1e-6, rtol=1e-6)
    @test opnorm(rho_one - exact_occupied_projector(H, 1)) < 2e-3

    shifted_sys, shifted_params, shifted_H = mcweeny_mu_system(shift=10.0)
    shifted_bounds = exact_spectral_bounds(shifted_H; padding=0.5)
    shifted = mcweeny_mu_result(shifted_sys, shifted_params, shifted_bounds, 10.0)
    @test shifted.converged
    @test opnorm(dense_matrix(shifted.rho, shifted_sys) - rho_half) < 2e-3

    @test_throws ArgumentError construct_rho_0(
        sys, params, bounds...; method=:mcweeny_mu,
    )
    @test_throws ArgumentError perform_purification(
        half_filled.rho, params; method=:mcweeny_mu, spectral_bounds=bounds,
    )
end
