@testset "P2.1 non-interacting square SCF" begin
    params = parameters_square(
        L=4,
        t=(-0.6, -0.35),
        W=(x, y) -> 0.11x + 0.07y + 0.013x * y,
        U=0.0,
        S=nothing,
        density=0.5,
        purification_steps=35,
        itensors_maxdim=64,
        scf_max_iterations=5,
    )
    for method in (:sp2, :palser_manolopoulos)
        sys = System(params)
        H0 = dense_matrix(sys.H0, sys)
        bounds = exact_spectral_bounds(H0; padding=0.5)
        exact = exact_occupied_projector(H0, 8)
        @test run_scf!(
            sys, bounds...;
            purification_method=method,
            verify_spectral_bounds=true,
            verbose=:nothing,
        )
        rho = dense_matrix(sys.ρ, sys)
        vh = dense_matrix(sys.VH, sys)
        vf = dense_matrix(sys.VF, sys)
        @test opnorm(rho - exact) < 3e-3
        @test norm(vh) == 0.0
        @test norm(vf) == 0.0
        @test isapprox(real(tr(rho)), 8.0; atol=3e-3)
        @test opnorm(H0 * rho - rho * H0) < 1e-7
    end
end
