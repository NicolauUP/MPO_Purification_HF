function dense_open_chain_reference(N, hopping, potential=nothing)
    H = zeros(Float64, N, N)
    for site in 0:(N - 1)
        H[site + 1, site + 1] = isnothing(potential) ? 0.0 : potential(site)
    end
    for bond in 0:(N - 2)
        amplitude = hopping isa Number ? hopping : hopping(bond)
        H[bond + 1, bond + 2] = amplitude
        H[bond + 2, bond + 1] = amplitude
    end
    return H
end

@testset "M1.2 one-dimensional static Hamiltonian" begin
    constant_hopping = -0.7
    functional_hopping(x) = -0.2 - 0.07x
    uniform_potential(x) = 0.35
    staggered_potential(x) = iseven(x) ? 0.4 : -0.4

    for L in 1:3
        N = 2^L

        constant_system = System(parameters_1d(L=L, t=constant_hopping))
        constant_dense = dense_matrix(constant_system.H0, constant_system)
        @test isapprox(
            constant_dense,
            dense_open_chain_reference(N, constant_hopping);
            atol=1e-12,
            rtol=1e-12,
        )
        if N > 2
            @test abs(constant_dense[1, N]) <= 1e-12
            @test abs(constant_dense[N, 1]) <= 1e-12
        end

        functional_system = System(parameters_1d(L=L, t=functional_hopping))
        functional_dense = dense_matrix(functional_system.H0, functional_system)
        @test isapprox(
            functional_dense,
            dense_open_chain_reference(N, functional_hopping);
            atol=1e-10,
            rtol=1e-10,
        )
        @test isapprox(functional_dense, functional_dense'; atol=1e-12, rtol=1e-12)

        uniform_system = System(parameters_1d(
            L=L, t=constant_hopping, W=uniform_potential,
        ))
        @test isapprox(
            dense_matrix(uniform_system.H0, uniform_system),
            dense_open_chain_reference(N, constant_hopping, uniform_potential);
            atol=1e-10,
            rtol=1e-10,
        )

        staggered_system = System(parameters_1d(
            L=L, t=constant_hopping, W=staggered_potential,
        ))
        @test isapprox(
            dense_matrix(staggered_system.H0, staggered_system),
            dense_open_chain_reference(N, constant_hopping, staggered_potential);
            atol=1e-10,
            rtol=1e-10,
        )
    end

    # The analytic cosine path must agree coefficient-for-coefficient with
    # the established generic function/TCI construction, while avoiding any
    # interpolation search. This also verifies the bond-index convention:
    # t(i) weights the open-chain bond i -> i+1.
    for L in 1:5
        N = 2^L
        cosine = CosineHopping(-1.0, -0.6, 0.73, 0.37)
        analytic_system = System(parameters_1d(L=L, t=cosine, U=0.0))
        tci_system = System(parameters_1d(L=L, t=x -> cosine(x), U=0.0))
        analytic_dense = dense_matrix(analytic_system.H0, analytic_system)
        tci_dense = dense_matrix(tci_system.H0, tci_system)
        reference = dense_open_chain_reference(N, cosine)

        @test isapprox(analytic_dense, reference; atol=2e-12, rtol=2e-12)
        @test isapprox(analytic_dense, tci_dense; atol=2e-10, rtol=2e-10)
        @test isapprox(analytic_dense, analytic_dense'; atol=1e-12, rtol=1e-12)
    end
end
