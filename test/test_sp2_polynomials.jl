@testset "M3.3 SP2 scalar polynomial invariants" begin
    for x in range(0.0, 1.0; length=101)
        square = x^2
        hole = 2x - x^2
        @test 0.0 <= square <= 1.0
        @test 0.0 <= hole <= 1.0
        @test square <= x
        @test hole >= x
    end

    @test 0.0^2 == 0.0
    @test 1.0^2 == 1.0
    @test 2 * 0.0 - 0.0^2 == 0.0
    @test 2 * 1.0 - 1.0^2 == 1.0

    params = parameters_1d(
        t=0.0,
        W=x -> (-2.0, -0.5, 0.7, 2.0)[Int(x) + 1],
        U=0.0,
        density=0.5,
    )
    sys = System(params)
    H = dense_matrix(sys.H0, sys)
    lambda_min, lambda_max = extrema(eigvals(Hermitian(H)))
    @test_throws ArgumentError construct_rho_0(
        sys, params, lambda_min + 0.1, lambda_max + 0.5;
        method=:sp2,
        verify_spectral_bounds=true,
    )
end
