@testset "P1.5 geometry-specific mean-field dispatch" begin
    sys = System(parameters_square(U=0.3))
    sys.ρ = Identity_MPO(sys.sites) * 0.5
    vh, vf = extract_mean_fields(sys)
    expected_vh = extract_hartree_mpo_square(sys)
    expected_vf = +(
        extract_fock_mpo_square_horizontal(sys),
        extract_fock_mpo_square_vertical(sys);
        cutoff=sys.params.itensors_tol,
        maxdim=sys.params.itensors_maxdim,
    )
    @test opnorm(dense_matrix(vh, sys) - dense_matrix(expected_vh, sys)) < 1e-12
    @test opnorm(dense_matrix(vf, sys) - dense_matrix(expected_vf, sys)) < 1e-12
end
