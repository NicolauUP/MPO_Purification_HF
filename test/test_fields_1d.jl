@testset "M4.1 1D Hartree/Fock dense reference" begin
    params = parameters_1d(L=2, U=0.7, S=nothing)
    sys = System(params)
    occupations = (0.15, 0.40, 0.65, 0.25)
    bonds = (0.03, -0.07, 0.11)
    diagonal = MPO_MeanField.diagonal_mpo_from_function(
        x -> occupations[Int(x) + 1], Float64, sys.sites, params.tci_tol,
    )
    bond_diagonal = MPO_MeanField.diagonal_mpo_from_function(
        x -> Int(x) < 3 ? bonds[Int(x) + 1] : 0.0, Float64, sys.sites, params.tci_tol,
    )
    T_R, T_L = MPO_MeanField.build_translation_chain(sys.sites)
    rho_right = apply(bond_diagonal, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    rho_left = apply(T_L, ITensors.dag(bond_diagonal); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    sys.ρ = +(diagonal, rho_right, rho_left; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

    rho = dense_matrix(sys.ρ, sys)
    cached_right, cached_left = sys.translations
    fresh_right, fresh_left = MPO_MeanField.build_translation_chain(sys.sites)
    @test dense_matrix(cached_right, sys) == dense_matrix(fresh_right, sys)
    @test dense_matrix(cached_left, sys) == dense_matrix(fresh_left, sys)
    vh = dense_matrix(extract_hartree_mpo_1d(sys), sys)
    vf = dense_matrix(extract_fock_mpo_1d(sys), sys)
    expected_vh = zeros(4, 4)
    expected_vf = zeros(4, 4)
    for i in 1:4
        for j in (i - 1, i + 1)
            1 <= j <= 4 && (expected_vh[i, i] += params.U * rho[j, j])
        end
        if i < 4
            expected_vf[i, i + 1] = -params.U * real(rho[i, i + 1])
            expected_vf[i + 1, i] = expected_vf[i, i + 1]
        end
    end
    @test opnorm(vh - expected_vh) < 1e-9
    @test opnorm(vf - expected_vf) < 1e-9
    @test opnorm(vh - vh') < 1e-12
    @test opnorm(vf - vf') < 1e-12

    evaluator = MPO_MeanField.HartreeEvaluate1D(sys, Dict{Int,Float64}())
    evaluator(0.0)
    evaluator(1.0)
    @test sort(collect(keys(evaluator.density_cache))) == [1, 2, 3]

    tensorial_hartree = dense_matrix(extract_hartree_mpo_tensorial_1d(sys), sys)
    @test tensorial_hartree ≈ vh atol=1e-10

    carry_hartree = dense_matrix(extract_hartree_mpo_binary_carry_1d(sys), sys)
    @test carry_hartree ≈ vh atol=1e-10
end
