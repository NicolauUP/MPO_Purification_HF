@testset "McWeeny standard and Horner polynomial forms" begin
    dense_rho = [
        0.25 0.03 0.00 0.01
        0.03 0.45 0.02 0.00
        0.00 0.02 0.65 0.04
        0.01 0.00 0.04 0.85
    ]
    dense_square = dense_rho * dense_rho
    dense_standard = 3dense_square - 2dense_square * dense_rho
    dense_horner = dense_square * (3I - 2dense_rho)
    @test dense_standard ≈ dense_horner atol=1e-14 rtol=1e-14

    params = parameters_1d(
        t=0.0,
        W=x -> (-1.5, -0.4, 0.6, 1.4)[Int(x) + 1],
        U=0.0,
        density=0.5,
        itensors_tol=1e-14,
        itensors_maxdim=64,
    )
    sys = System(params)
    rho = construct_rho_0(
        sys, params, -2.0, 2.0;
        method=:mcweeny_mu,
        chemical_potential=0.0,
    )
    rho_squared = apply(rho, rho;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    standard = MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:standard,
    )
    horner = MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params;
        form=:horner,
        identity_mpo=Identity_MPO(sys.sites),
    )

    @test standard.auxiliary_kind == :rho_cubed
    @test horner.auxiliary_kind == :horner_factor
    @test opnorm(dense_matrix(standard.rho, sys) - dense_matrix(horner.rho, sys)) < 1e-11
    @test_throws ArgumentError MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:horner,
    )
    @test_throws ArgumentError MPO_MeanField._mcweeny_update_from_square(
        rho, rho_squared, params; form=:unknown,
    )
end
