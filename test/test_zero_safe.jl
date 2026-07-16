@testset "M0 zero-safe MPO paths" begin
    sys_none = System(parameters_1d(S=nothing))
    @test norm(dense_matrix(sys_none.VH, sys_none)) == 0.0

    sys_zero = System(parameters_1d(S=x -> 0.0))
    @test norm(dense_matrix(sys_zero.VH, sys_zero)) == 0.0

    square_none = System(parameters_square(S=nothing))
    @test norm(dense_matrix(square_none.VH, square_none)) == 0.0

    noninteracting = System(parameters_1d(U=0.0, S=nothing))
    noninteracting.ρ = Identity_MPO(noninteracting.sites) * 0.5
    vh = extract_hartree_mpo_1d(noninteracting)
    vf = extract_fock_mpo_1d(noninteracting)
    @test norm(dense_matrix(vh, noninteracting)) == 0.0
    @test norm(dense_matrix(vf, noninteracting)) == 0.0

    zero_density = System(parameters_1d(U=0.3, S=nothing))
    vh_zero_density = extract_hartree_mpo_1d(zero_density)
    vf_zero_density = extract_fock_mpo_1d(zero_density)
    @test norm(dense_matrix(vh_zero_density, zero_density)) == 0.0
    @test norm(dense_matrix(vf_zero_density, zero_density)) == 0.0

    scf_system = System(parameters_1d(
        U=0.0,
        S=nothing,
        purification_steps=20,
        scf_max_iterations=3,
    ))
    @test run_scf!(scf_system, -3.0, 3.0; verbose=:nothing)
    @test norm(dense_matrix(scf_system.VH, scf_system)) == 0.0
    @test norm(dense_matrix(scf_system.VF, scf_system)) == 0.0

    @test MPO_MeanField.safe_relative_change(0.0, 0.0) == 0.0
    @test isfinite(MPO_MeanField.safe_relative_change(1.0, 0.0))
    @test_throws ArgumentError MPO_MeanField.safe_relative_change(-1.0, 1.0)
end
