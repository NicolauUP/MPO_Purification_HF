function _set_real_tridiagonal_density!(sys, diagonal, bonds)
    params = sys.params
    diagonal_mpo = MPO_MeanField.diagonal_mpo_from_function(
        x -> diagonal[Int(x) + 1], Float64, sys.sites, params.tci_tol,
    )
    bond_mpo = MPO_MeanField.diagonal_mpo_from_function(
        x -> Int(x) < length(bonds) ? bonds[Int(x) + 1] : 0.0,
        Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L = MPO_MeanField.build_translation_chain(sys.sites)
    rho_right = apply(bond_mpo, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_left = apply(T_L, ITensors.dag(bond_mpo); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.rho = +(diagonal_mpo, rho_right, rho_left; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    return sys
end

@testset "M6.1 nearest-neighbour 1D HF energy" begin
    params = parameters_1d(L=2, t=-0.6, U=0.7, W=x -> (0.2, -0.1, 0.3, -0.2)[Int(x) + 1])
    sys = System(params)
    diagonal = [0.15, 0.40, 0.65, 0.25]
    bonds = [0.03, -0.07, 0.11]
    _set_real_tridiagonal_density!(sys, diagonal, bonds)

    rho = dense_matrix(sys.rho, sys)
    h0 = dense_matrix(sys.H0, sys)
    energy = nearest_neighbor_hf_energy_1d(sys)
    expected_hartree = params.U * sum(diagonal[i] * diagonal[i + 1] for i in 1:3)
    expected_fock = -params.U * sum(abs2(bonds[i]) for i in 1:3)
    @test isapprox(energy.kinetic, real(tr(h0 * rho)); atol=1e-10)
    @test isapprox(energy.hartree, expected_hartree; atol=1e-10)
    @test isapprox(energy.fock, expected_fock; atol=1e-10)
    @test isapprox(energy.total, energy.kinetic + energy.hartree + energy.fock; atol=1e-10)

    vh = dense_matrix(extract_hartree_mpo_1d(sys), sys)
    vf = dense_matrix(extract_fock_mpo_1d(sys), sys)
    @test isapprox(energy.interaction, real(tr((vh + vf) * rho)) / 2; atol=1e-10)

    # For real, symmetric density variations, the finite difference of the
    # implemented interaction functional is the extracted mean field.
    direction_diagonal = [0.12, -0.08, 0.05, -0.09]
    direction_bonds = [0.04, -0.02, 0.03]
    direction = zeros(4, 4)
    for i in 1:4
        direction[i, i] = direction_diagonal[i]
    end
    for i in 1:3
        direction[i, i + 1] = direction_bonds[i]
        direction[i + 1, i] = direction_bonds[i]
    end
    epsilon = 1e-6
    _set_real_tridiagonal_density!(sys, diagonal .+ epsilon .* direction_diagonal, bonds .+ epsilon .* direction_bonds)
    e_plus = nearest_neighbor_hf_energy_1d(sys).interaction
    _set_real_tridiagonal_density!(sys, diagonal .- epsilon .* direction_diagonal, bonds .- epsilon .* direction_bonds)
    e_minus = nearest_neighbor_hf_energy_1d(sys).interaction
    finite_difference = (e_plus - e_minus) / (2epsilon)

    _set_real_tridiagonal_density!(sys, diagonal, bonds)
    vh = dense_matrix(extract_hartree_mpo_1d(sys), sys)
    vf = dense_matrix(extract_fock_mpo_1d(sys), sys)
    @test isapprox(finite_difference, real(tr((vh + vf) * direction)); atol=1e-8)
end
