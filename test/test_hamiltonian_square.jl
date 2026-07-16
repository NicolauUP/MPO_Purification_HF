function square_site_index(x, y, L)
    index = 0
    for bit in 0:(div(L, 2) - 1)
        index |= ((y >> bit) & 1) << (2bit)
        index |= ((x >> bit) & 1) << (2bit + 1)
    end
    return index + 1
end

function dense_square_translation_reference(L, direction)
    side = 2^div(L, 2)
    T = zeros(Float64, side^2, side^2)
    for x in 0:(side - 1), y in 0:(side - 1)
        row = square_site_index(x, y, L)
        if direction == :right && x < side - 1
            T[row, square_site_index(x + 1, y, L)] = 1.0
        elseif direction == :left && x > 0
            T[row, square_site_index(x - 1, y, L)] = 1.0
        elseif direction == :up && y < side - 1
            T[row, square_site_index(x, y + 1, L)] = 1.0
        elseif direction == :down && y > 0
            T[row, square_site_index(x, y - 1, L)] = 1.0
        end
    end
    return T
end

function dense_square_hamiltonian_reference(L, tx, ty, potential=nothing)
    side = 2^div(L, 2)
    H = zeros(Float64, side^2, side^2)
    hopping_x(x, y) = tx isa Number ? tx : tx(x, y)
    hopping_y(x, y) = ty isa Number ? ty : ty(x, y)
    for x in 0:(side - 1), y in 0:(side - 1)
        site = square_site_index(x, y, L)
        H[site, site] = isnothing(potential) ? 0.0 : potential(x, y)
        if x < side - 1
            neighbour = square_site_index(x + 1, y, L)
            H[site, neighbour] = H[neighbour, site] = hopping_x(x, y)
        end
        if y < side - 1
            neighbour = square_site_index(x, y + 1, L)
            H[site, neighbour] = H[neighbour, site] = hopping_y(x, y)
        end
    end
    return H
end

@testset "M1.3 square static Hamiltonian" begin
    L = 4
    sites = siteinds("Qubit", L)
    bra_states, ket_states = MPO_MeanField.precompute_qtt_states(sites)
    dense_translation(mpo) = [
        MatrixChecker(mpo, sites, i, j, bra_states, ket_states)
        for i in 1:16, j in 1:16
    ]
    T_R, T_L, T_U, T_D = build_translation_square(sites)

    @test dense_translation(T_R) == dense_square_translation_reference(L, :right)
    @test dense_translation(T_L) == dense_square_translation_reference(L, :left)
    @test dense_translation(T_U) == dense_square_translation_reference(L, :up)
    @test dense_translation(T_D) == dense_square_translation_reference(L, :down)
    @test dense_translation(T_L) == dense_translation(T_R)'
    @test dense_translation(T_D) == dense_translation(T_U)'

    constant = System(parameters_square(L=L, t=(-0.7, -0.3)))
    @test isapprox(
        dense_matrix(constant.H0, constant),
        dense_square_hamiltonian_reference(L, -0.7, -0.3);
        atol=1e-12, rtol=1e-12,
    )

    tx(x, y) = -0.2 - 0.03x + 0.01y
    ty(x, y) = -0.4 + 0.02x - 0.05y
    W(x, y) = 0.1x - 0.2y
    functional = System(parameters_square(L=L, t=(tx, ty), W=W))
    @test isapprox(
        dense_matrix(functional.H0, functional),
        dense_square_hamiltonian_reference(L, tx, ty, W);
        atol=1e-10, rtol=1e-10,
    )
    @test isapprox(
        dense_matrix(functional.H0, functional),
        dense_matrix(functional.H0, functional)';
        atol=1e-10, rtol=1e-10,
    )

    mixed = System(parameters_square(L=L, t=(-0.6, ty), W=W))
    @test isapprox(
        dense_matrix(mixed.H0, mixed),
        dense_square_hamiltonian_reference(L, -0.6, ty, W);
        atol=1e-10, rtol=1e-10,
    )
end
