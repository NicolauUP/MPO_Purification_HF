@testset "square interacting SCF spectral bounds" begin
    constant = parameters_square(
        L=2, t=(-0.6, -0.35), U=0.15, W=nothing, S=nothing,
        purification_steps=35, itensors_maxdim=64, scf_mixing=0.4,
        scf_tol=0.1, scf_max_iterations=40,
    )
    constant_bounds = square_scf_spectral_bounds(constant; margin=0.5)
    @test isapprox(constant_bounds[1], -3.6; atol=1e-12)
    @test isapprox(constant_bounds[2], 3.6; atol=1e-12)

    potential(x, y) = 0.11x + 0.07y + 0.013x * y
    params = parameters_square(
        L=2, t=(-0.6, -0.35), U=0.15, W=potential, S=nothing,
        purification_steps=35, itensors_maxdim=64, scf_mixing=0.4,
        scf_tol=0.1, scf_max_iterations=40,
    )
    bounds = square_scf_spectral_bounds(
        params;
        potential_bounds=(0.0, potential(1, 1)),
        margin=0.5,
    )
    @test isapprox(bounds[1], -3.6; atol=1e-12)
    @test isapprox(bounds[2], 3.793; atol=1e-12)
    @test_throws ArgumentError square_scf_spectral_bounds(params)
    @test_throws ArgumentError square_scf_spectral_bounds(
        params; potential_bounds=(1.0, -1.0),
    )

    functional_hopping = parameters_square(
        L=2, t=((x, y) -> -0.6, -0.35), U=0.15, W=nothing, S=nothing,
        purification_steps=20, itensors_maxdim=64, scf_max_iterations=5,
    )
    @test_throws ArgumentError square_scf_spectral_bounds(functional_hopping)
    functional_bounds = square_scf_spectral_bounds(
        functional_hopping; hopping_abs_bounds=(0.6, 0.35), margin=0.5,
    )
    @test isapprox(functional_bounds[1], -3.6; atol=1e-12)
    @test isapprox(functional_bounds[2], 3.6; atol=1e-12)

    sys = System(params)
    @test run_scf!(
        sys, bounds...;
        purification_method=:sp2,
        verify_spectral_bounds=true,
        verbose=:nothing,
    )
    history = scf_diagnostics(sys).history
    @test all(record -> !isnothing(record.rho_bond_dimension), history)
    @test all(record -> !isnothing(record.hartree_bond_dimension), history)
    @test all(record -> !isnothing(record.fock_bond_dimension), history)
    @test all(record -> !isnothing(record.effective_hamiltonian_bond_dimension), history)
    @test all(record -> record.rho_bond_dimension <= params.itensors_maxdim, history)
    H_effective = +(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    extrema = verify_spectral_bounds_exact(sys, H_effective, bounds...)
    @test bounds[1] <= extrema[1] <= extrema[2] <= bounds[2]
end
