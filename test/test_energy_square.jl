function _set_real_square_density!(sys, diagonal, horizontal, vertical)
    params = sys.params
    diagonal_mpo = MPO_MeanField.diagonal_mpo_from_function(
        coordinate -> diagonal[Int(coordinate) + 1], Float64, sys.sites, params.tci_tol,
    )
    horizontal_mpo = MPO_MeanField.diagonal_mpo_from_function(
        coordinate -> begin
            x, y = square_lattice_decoder(Int(coordinate), params.L)
            x < 2^div(params.L, 2) - 1 ? horizontal[square_lattice_index(x, y, params.L)] : 0.0
        end,
        Float64,
        sys.sites,
        params.tci_tol,
    )
    vertical_mpo = MPO_MeanField.diagonal_mpo_from_function(
        coordinate -> begin
            x, y = square_lattice_decoder(Int(coordinate), params.L)
            y < 2^div(params.L, 2) - 1 ? vertical[square_lattice_index(x, y, params.L)] : 0.0
        end,
        Float64,
        sys.sites,
        params.tci_tol,
    )
    T_R, T_L, T_U, T_D = sys.translations
    rho_right = apply(horizontal_mpo, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_left = apply(T_L, ITensors.dag(horizontal_mpo); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_up = apply(vertical_mpo, T_U; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_down = apply(T_D, ITensors.dag(vertical_mpo); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.ρ = +(diagonal_mpo, rho_right, rho_left, rho_up, rho_down;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    return sys
end

@testset "P3.2 nearest-neighbour square HF energy" begin
    params = parameters_square(L=4, t=(-0.6, -0.35), U=0.7,
        W=(x, y) -> 0.1x - 0.07y, S=nothing)
    sys = System(params)
    diagonal = [0.12 + 0.015 * site for site in 1:16]
    horizontal = [0.01 * (1 + site) for site in 1:16]
    vertical = [-0.008 * (1 + site) for site in 1:16]
    _set_real_square_density!(sys, diagonal, horizontal, vertical)

    rho = dense_matrix(sys.ρ, sys)
    h0 = dense_matrix(sys.H0, sys)
    energy = nearest_neighbor_hf_energy_square(sys)
    energy_with_explicit_h0 = nearest_neighbor_hf_energy_square(
        sys; hopping_hamiltonian=sys.H0,
    )
    expected_hartree = 0.0
    expected_fock = 0.0
    for (site, neighbour, _) in square_undirected_bonds(params.L)
        expected_hartree += params.U * diagonal[site] * diagonal[neighbour]
        expected_fock -= params.U * real(rho[site, neighbour])^2
    end
    @test isapprox(energy.kinetic, real(tr(h0 * rho)); atol=1e-10)
    @test isapprox(energy.hartree, expected_hartree; atol=1e-10)
    @test isapprox(energy.fock, expected_fock; atol=1e-10)
    @test isapprox(energy.interaction, energy.hartree + energy.fock; atol=1e-10)
    @test isapprox(energy.total, energy.kinetic + energy.interaction; atol=1e-10)
    @test energy_with_explicit_h0 == energy

    vh, vf = extract_mean_fields(sys)
    mean_field = dense_matrix(+(vh, vf;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    ), sys)
    @test isapprox(energy.interaction, real(tr(mean_field * rho)) / 2; atol=1e-10)
end

@testset "P3.3 square observables" begin
    params = parameters_square(L=4, t=(-0.6, -0.35), U=0.7, W=nothing, S=nothing)
    sys = System(params)
    diagonal = [0.12 + 0.015 * site for site in 1:16]
    horizontal = [0.01 * (1 + site) for site in 1:16]
    vertical = [-0.008 * (1 + site) for site in 1:16]
    _set_real_square_density!(sys, diagonal, horizontal, vertical)
    sys.VH, sys.VF = extract_mean_fields(sys)

    result = observables_square(sys)
    generic_result = observables(sys)
    rho = dense_matrix(sys.ρ, sys)
    h_effective = dense_matrix(+(sys.H0, sys.VH, sys.VF;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    ), sys)
    expected_horizontal = Tuple{Int,Int}[]
    expected_vertical = Tuple{Int,Int}[]
    for (site, neighbour, orientation) in square_undirected_bonds(params.L)
        if orientation == :horizontal
            push!(expected_horizontal, (site, neighbour))
        else
            @test orientation == :vertical
            push!(expected_vertical, (site, neighbour))
        end
    end

    @test result.site_density ≈ diagonal atol=1e-10
    @test result.horizontal_bonds == expected_horizontal
    @test result.vertical_bonds == expected_vertical
    @test result.horizontal_bond_order ≈ [rho[site, neighbour] for (site, neighbour) in expected_horizontal] atol=1e-10
    @test result.vertical_bond_order ≈ [rho[site, neighbour] for (site, neighbour) in expected_vertical] atol=1e-10
    @test isapprox(result.particle_number, sum(diagonal); atol=1e-10)
    @test result.energy == nearest_neighbor_hf_energy_square(sys)
    @test generic_result == result
    @test result.hermiticity_residual < 1e-10
    @test isapprox(result.idempotency_residual, norm(rho * rho - rho) / norm(rho); atol=1e-10)
    @test isapprox(result.stationarity_residual, norm(h_effective * rho - rho * h_effective) / norm(h_effective * rho); atol=1e-10)
end
