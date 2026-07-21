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

    carry_hartree = dense_matrix(extract_hartree_mpo_binary_carry_square(sys), sys)
    @test carry_hartree ≈ expected atol=1e-10
    @test opnorm(carry_hartree - carry_hartree') < 1e-12

    adjacency_hartree = dense_matrix(
        extract_hartree_mpo_binary_carry_square_adjacency(sys), sys,
    )
    @test adjacency_hartree ≈ expected atol=1e-10
    @test opnorm(adjacency_hartree - adjacency_hartree') < 1e-12

    # The square SCF selector must use the fused adjacency path, not the
    # legacy TCI Hartree extractor. Check every site directly against the
    # physical nearest-neighbour density sum.
    scf_hartree, _ = extract_mean_fields(sys)
    for site in 1:16
        direct = params.U * sum(
            real(MatrixChecker(
                sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states,
            ))
            for neighbour in values(square_neighbours(site, params.L))
            if !isnothing(neighbour)
        )
        observed = real(MatrixChecker(
            scf_hartree, sys.sites, site, site, sys.bra_states, sys.ket_states,
        ))
        @test isapprox(observed, direct; atol=1e-10, rtol=1e-10)
    end

    # Check each interleaved binary-carry direction separately. The final
    # Hartree sum alone could conceal a horizontal/vertical axis swap.
    density_tensors = MPO_MeanField._density_diagonal_qtt_tensors(sys.ρ, sys.sites)
    carry_components = (
        right=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :right),
            sys.sites, params; symmetrize=false,
        ),
        left=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :left),
            sys.sites, params; symmetrize=false,
        ),
        up=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :up),
            sys.sites, params; symmetrize=false,
        ),
        down=MPO_MeanField._diagonal_mpo_from_qtt_tensors(
            MPO_MeanField._shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :down),
            sys.sites, params; symmetrize=false,
        ),
    )
    for site in 1:16
        neighbours = square_neighbours(site, params.L)
        for direction in (:right, :left, :up, :down)
            neighbour = getproperty(neighbours, direction)
            expected_density = isnothing(neighbour) ? 0.0 : occupations[neighbour]
            observed_density = MatrixChecker(
                getproperty(carry_components, direction), sys.sites, site, site,
                sys.bra_states, sys.ket_states,
            )
            @test isapprox(observed_density, expected_density; atol=1e-10)
        end
    end
end

@testset "P1.2b square adjacency Hartree with exact checkerboard density" begin
    params = parameters_square(L=4, U=0.7, S=nothing, itensors_tol=1e-12)
    sys = System(params)
    delta = 0.2
    density = (x, y) -> 0.5 + (iseven(x + y) ? delta : -delta)
    sys.ρ = MPO_MeanField.diagonal_mpo_from_function(
        z -> begin
            x, y = square_lattice_decoder(Int(z), params.L)
            density(x, y)
        end,
        Float64,
        sys.sites,
        params.tci_tol,
    )

    expected = zeros(16, 16)
    for site in 1:16
        expected[site, site] = params.U * sum(
            density(square_lattice_decoder(neighbour - 1, params.L)...)
            for neighbour in values(square_neighbours(site, params.L))
            if !isnothing(neighbour)
        )
    end
    four_carry = dense_matrix(extract_hartree_mpo_binary_carry_square(sys), sys)
    adjacency = dense_matrix(extract_hartree_mpo_binary_carry_square_adjacency(sys), sys)
    @test four_carry ≈ expected atol=1e-10
    @test adjacency ≈ expected atol=1e-10
    @test adjacency ≈ four_carry atol=1e-10
end

@testset "P1.2a square binary-carry Hartree with checkerboard SP2 density" begin
    # Static checkerboard CDW potential, not an SCF seed or calculation.
    W = (x, y) -> iseven(Int(x) + Int(y)) ? 0.6 : -0.6
    params = parameters_square(
        # Keep the SP2-backed unit test small. At larger L the current generic
        # four-MPO addition is itself the performance subject under study.
        L=2,
        t=(-0.6, -0.35),
        U=0.7,
        W=W,
        S=nothing,
        itensors_maxdim=128,
        purification_steps=60,
    )
    sys = System(params)
    spectral_bounds = (-3.5, 3.5)
    rho0 = construct_rho_0(sys, params, spectral_bounds...; method=:sp2)
    result = perform_purification(
        rho0,
        params;
        method=:sp2,
        verbose=0,
        spectral_bounds=spectral_bounds,
    )
    @test result.converged
    sys.ρ = result.rho

    carry_hartree = extract_hartree_mpo_binary_carry_square(sys)
    for site in 1:(2^params.L)
        direct_hartree = params.U * sum(
            real(MatrixChecker(
                sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states,
            ))
            for neighbour in values(square_neighbours(site, params.L))
            if !isnothing(neighbour)
        )
        carry_value = real(MatrixChecker(
            carry_hartree, sys.sites, site, site, sys.bra_states, sys.ket_states,
        ))
        @test isapprox(carry_value, direct_hartree; atol=1e-8, rtol=1e-8)
    end
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
