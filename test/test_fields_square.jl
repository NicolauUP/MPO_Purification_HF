@testset "P1.2 square Hartree dense reference" begin
    params = parameters_square(L=4, U=0.7, S=nothing)
    sys = System(params)
    occupations = [0.05 + 0.03 * i for i in 1:16]
    sys.ρ = MPO_MeanField.diagonal_mpo_from_function(
        x -> occupations[Int(x) + 1], Float64, sys.sites, params.tci_tol,
    )

    hartree = dense_matrix(extract_hartree_mpo_square(sys), sys)
    expected = zeros(16, 16)
    for site in 1:16
        for neighbour in values(square_neighbours(site, params.L))
            !isnothing(neighbour) && (expected[site, site] += params.U * occupations[neighbour])
        end
    end
    @test opnorm(hartree - expected) < 1e-10
    @test opnorm(hartree - hartree') < 1e-12

    evaluator = MPO_MeanField.HartreeEvaluateSquare(sys, Dict{Int,Float64}())
    evaluator(0.0)
    evaluator(1.0)
    expected_cached_sites = Int[
        neighbour for site in (1, 2)
        for neighbour in values(square_neighbours(site, params.L))
        if !isnothing(neighbour)
    ]
    @test sort(collect(keys(evaluator.density_cache))) == sort(unique(expected_cached_sites))

    tensorial_hartree = dense_matrix(extract_hartree_mpo_tensorial_square(sys), sys)
    @test tensorial_hartree ≈ hartree atol=1e-10
end

@testset "P1.3 square horizontal Fock dense reference" begin
    params = parameters_square(L=4, U=0.7, S=nothing)
    sys = System(params)
    occupations = x -> 0.1 + 0.01 * Int(x)
    horizontal_bond = coordinate -> begin
        x, y = square_lattice_decoder(Int(coordinate), params.L)
        return x < 3 ? 0.02 * (1 + x + 2y) : 0.0
    end
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        occupations, Float64, sys.sites, params.tci_tol,
    )
    bond_diagonal = MPO_MeanField.diagonal_mpo_from_function(
        horizontal_bond, Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L, _, _ = sys.translations
    rho_right = apply(bond_diagonal, T_R;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    rho_left = apply(T_L, ITensors.dag(bond_diagonal);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    sys.ρ = +(diagonal, rho_right, rho_left;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )

    rho = dense_matrix(sys.ρ, sys)
    fock = dense_matrix(extract_fock_mpo_square_horizontal(sys), sys)
    expected = zeros(16, 16)
    for (site, right, orientation) in square_undirected_bonds(params.L)
        orientation == :horizontal || continue
        expected[site, right] = -params.U * real(rho[site, right])
        expected[right, site] = expected[site, right]
    end
    @test opnorm(fock - expected) < 1e-10
    @test opnorm(fock - fock') < 1e-12
end

@testset "P1.4 square vertical Fock dense reference" begin
    params = parameters_square(L=4, U=0.7, S=nothing)
    sys = System(params)
    occupations = x -> 0.1 + 0.01 * Int(x)
    vertical_bond = coordinate -> begin
        x, y = square_lattice_decoder(Int(coordinate), params.L)
        return y < 3 ? 0.015 * (1 + 2x + y) : 0.0
    end
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        occupations, Float64, sys.sites, params.tci_tol,
    )
    bond_diagonal = MPO_MeanField.diagonal_mpo_from_function(
        vertical_bond, Float64, sys.sites, params.tci_tol,
    )
    _, _, T_U, T_D = sys.translations
    rho_up = apply(bond_diagonal, T_U;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    rho_down = apply(T_D, ITensors.dag(bond_diagonal);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    sys.ρ = +(diagonal, rho_up, rho_down;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )

    rho = dense_matrix(sys.ρ, sys)
    fock = dense_matrix(extract_fock_mpo_square_vertical(sys), sys)
    expected = zeros(16, 16)
    for (site, up, orientation) in square_undirected_bonds(params.L)
        orientation == :vertical || continue
        expected[site, up] = -params.U * real(rho[site, up])
        expected[up, site] = expected[site, up]
    end
    @test opnorm(fock - expected) < 1e-10
    @test opnorm(fock - fock') < 1e-12
end
