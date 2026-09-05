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

    carry_fock = dense_matrix(extract_fock_mpo_binary_carry_1d(sys), sys)
    @test carry_fock ≈ vf atol=1e-10
    @test opnorm(carry_fock - carry_fock') < 1e-12

    # The present mean-field convention deliberately retains only the real
    # exchange coefficient. This complex Hermitian density distinguishes a
    # direct superdiagonal contraction from a transposed/conjugated one.
    complex_bond_diagonal = (1 + 2im) * bond_diagonal
    complex_right = apply(complex_bond_diagonal, T_R;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    complex_left = apply(T_L, ITensors.dag(complex_bond_diagonal);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    sys.ρ = +(diagonal, complex_right, complex_left;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    complex_carry_fock = dense_matrix(extract_fock_mpo_binary_carry_1d(sys), sys)
    @test opnorm(complex_carry_fock - expected_vf) < 1e-10
    @test opnorm(complex_carry_fock - complex_carry_fock') < 1e-12
end

@testset "M4.2 1D binary-carry fields with SP2 projector" begin
    # Validate the public 1D field route on a genuinely purified density, not
    # only on a hand-built MPO. The staggered potential supplies a small gap.
    params = parameters_1d(
        L=2,
        t=-1.0,
        U=0.7,
        W=x -> iseven(Int(x)) ? 0.8 : -0.8,
        S=nothing,
        itensors_maxdim=64,
        purification_steps=40,
    )
    sys = System(params)
    bounds = (-3.0, 3.0)
    purified = perform_purification(
        construct_rho_0(sys, params, bounds...; method=:sp2),
        params;
        method=:sp2,
        verbose=0,
        spectral_bounds=bounds,
        sp2_idempotency_tolerance=2e-4,
        sp2_trace_tolerance=4e-6,
    )
    @test purified.converged
    sys.ρ = purified.rho

    legacy_hartree = extract_hartree_mpo_1d(sys)
    legacy_fock = extract_fock_mpo_1d(sys)
    carry_hartree, carry_fock = extract_mean_fields(sys)
    N = 1 << params.L
    for site in 1:N
        direct_hartree = params.U * sum(
            real(MatrixChecker(
                sys.ρ, sys.sites, neighbour, neighbour, sys.bra_states, sys.ket_states,
            ))
            for neighbour in (site - 1, site + 1)
            if 1 <= neighbour <= N
        )
        for field in (legacy_hartree, carry_hartree)
            observed = real(MatrixChecker(
                field, sys.sites, site, site, sys.bra_states, sys.ket_states,
            ))
            @test isapprox(observed, direct_hartree; atol=1e-10, rtol=1e-10)
        end
    end
    for site in 1:(N - 1)
        direct_fock = -params.U * real(MatrixChecker(
            sys.ρ, sys.sites, site, site + 1, sys.bra_states, sys.ket_states,
        ))
        for field in (legacy_fock, carry_fock)
            observed = real(MatrixChecker(
                field, sys.sites, site, site + 1, sys.bra_states, sys.ket_states,
            ))
            @test isapprox(observed, direct_fock; atol=1e-10, rtol=1e-10)
        end
    end
end
