@testset "QTT Hadamard MPS KPM probe" begin
    levels = 4
    sites = siteinds("Qubit", levels)
    row = 5
    probe = MPO_MeanField._qtt_hadamard_probe_mps(sites, row)
    expected_probe = MPO_MeanField._qtt_hadamard_probe_vector(levels, row)
    @test MPO_MeanField._qtt_mps_amplitudes(probe, sites) == expected_probe
    @test maxlinkdim(probe) == 1

    # A diagonal MPO gives a small deterministic check of recurrence
    # coefficients, amplitudes, and gauge, independent of the lattice MPO.
    diagonal = collect(range(-0.9, 0.9; length=1 << levels))
    _, H, _ = MPO_MeanField.Quantics_TCI(
        index -> diagonal[Int(index) + 1], Float64, sites, 1e-14,
    )
    coefficients = MPO_MeanField._kpm_coefficients(32, 0.0)
    mps_result = MPO_MeanField._qtt_mps_chebyshev_apply(
        H, probe, coefficients; cutoff=1e-12, maxdim=16,
    )
    function scalar_kpm(value)
        previous, current = 1.0, value
        result = coefficients[1] * previous / 2 + coefficients[2] * current
        for order in 2:(length(coefficients) - 1)
            following = 2 * value * current - previous
            result += coefficients[order + 1] * following
            previous, current = current, following
        end
        result
    end
    expected = expected_probe .* scalar_kpm.(diagonal)
    @test isapprox(
        MPO_MeanField._qtt_mps_amplitudes(mps_result.state, sites), expected;
        atol=2e-10,
        rtol=2e-10,
    )
    @test maximum(record.state_max_chi for record in mps_result.trajectory) <= 16

    # The finite-probe density estimator is a direct MPS sum. No TCI sampling
    # is needed because multiplying by a Walsh--Hadamard row is site-local.
    rows = [0, row]
    states = MPS[]
    for hadamard_row in rows
        row_probe = MPO_MeanField._qtt_hadamard_probe_mps(sites, hadamard_row)
        push!(states, MPO_MeanField._qtt_mps_chebyshev_apply(
            H, row_probe, coefficients; cutoff=1e-12, maxdim=16,
        ).state)
    end
    direct_density = MPO_MeanField._qtt_density_mps_from_hadamard_probes(
        states, sites, rows; cutoff=1e-12, maxdim=16,
    )
    raw_density = [sum(
        MPO_MeanField._qtt_hadamard_probe_vector(levels, hadamard_row)[index + 1] *
        MPO_MeanField._qtt_mps_amplitude(state, sites, index)
        for (hadamard_row, state) in zip(rows, states)
    ) / length(rows) for index in 0:((1 << levels) - 1)]
    @test isapprox(MPO_MeanField._qtt_mps_amplitudes(direct_density, sites), raw_density;
        atol=2e-10, rtol=2e-10)
    @test maxlinkdim(direct_density) <= 16

    # A block MPS carries a spectator probe index. MPO--MPS KPM propagation
    # must reproduce each independently propagated Hadamard column.
    block_rows = [0, row]
    block, probe_index = MPO_MeanField._qtt_hadamard_probe_block_mps(
        sites, block_rows; cutoff=1e-14, maxdim=128,
    )
    block_result = MPO_MeanField._qtt_mps_chebyshev_apply(
        H, block, coefficients; cutoff=1e-14, maxdim=128,
    )
    for (slot, state) in enumerate(states)
        expected_column = MPO_MeanField._qtt_mps_amplitudes(state, sites)
        actual_column = [MPO_MeanField._qtt_block_mps_amplitude(
            block_result.state, sites, probe_index, slot, index,
        ) for index in 0:((1 << levels) - 1)]
        @test isapprox(actual_column, expected_column; atol=3e-10, rtol=3e-10)
    end

    # A power-of-two Hadamard family can instead be represented by binary
    # probe bits. The first 2^q Gray rows are a permutation of the standard
    # Hadamard codes, so this changes only the column order, not the ensemble.
    register, register_sites, _ = MPO_MeanField._qtt_hadamard_probe_register_mps(sites, 2)
    @test maxlinkdim(register) <= 2
    for index in 0:((1 << levels) - 1), code in 0:3
        expected_amplitude = isodd(count_ones(index & code)) ? -1.0 : 1.0
        @test isapprox(MPO_MeanField._qtt_probe_register_amplitude(
            register, sites, register_sites, index, code,
        ), expected_amplitude; atol=1e-12, rtol=1e-12)
    end
    extended_H = MPO_MeanField._qtt_extend_mpo_with_probe_register(H, sites, register_sites)
    extended_result = MPO_MeanField._qtt_mps_chebyshev_apply(
        extended_H, register, coefficients; cutoff=1e-12, maxdim=32,
    )
    for index in 0:((1 << levels) - 1), code in 0:3
        expected_amplitude = (isodd(count_ones(index & code)) ? -1.0 : 1.0) *
            scalar_kpm(diagonal[index + 1])
        @test isapprox(MPO_MeanField._qtt_probe_register_amplitude(
            extended_result.state, sites, register_sites, index, code,
        ), expected_amplitude; atol=3e-10, rtol=3e-10)
    end

    # Contracting the register after the Walsh weighting gives the direct
    # density estimator as a spatial MPS, with no coordinate sampling.
    density_register = MPO_MeanField._qtt_density_mps_from_probe_register(
        extended_result.state, register, sites, register_sites;
        cutoff=1e-12, maxdim=32,
    )
    register_density = MPO_MeanField._qtt_mps_amplitudes(density_register, sites)
    exact_density = [sum(
        (isodd(count_ones(index & code)) ? -1.0 : 1.0) *
        MPO_MeanField._qtt_probe_register_amplitude(
            extended_result.state, sites, register_sites, index, code,
        ) for code in 0:3
    ) / 4 for index in 0:((1 << levels) - 1)]
    @test isapprox(register_density, exact_density; atol=3e-10, rtol=3e-10)

    # The existing binary translation MPO, lifted through the spectator
    # register, supplies the shifted Walsh column for a directed bond. This
    # checks the full local Fock-estimator contraction including open edges.
    T_R, _, _, _ = MPO_MeanField.build_translation_square(sites)
    right_bond_register = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        extended_result.state, register, T_R, sites, register_sites;
        cutoff=1e-12, maxdim=32,
    )
    register_right_bond = MPO_MeanField._qtt_mps_amplitudes(right_bond_register, sites)
    exact_right_bond = [begin
        x, y = MPO_MeanField.square_lattice_decoder(index, levels)
        if x == (1 << (levels ÷ 2)) - 1
            0.0
        else
            neighbour = MPO_MeanField.square_lattice_index(x + 1, y, levels) - 1
            sum(
                MPO_MeanField._qtt_probe_register_amplitude(
                    extended_result.state, sites, register_sites, index, code,
                ) * (isodd(count_ones(neighbour & code)) ? -1.0 : 1.0)
                for code in 0:3
            ) / 4
        end
    end for index in 0:((1 << levels) - 1)]
    @test isapprox(register_right_bond, exact_right_bond; atol=3e-10, rtol=3e-10)

    _, _, T_U, _ = MPO_MeanField.build_translation_square(sites)
    up_bond_register = MPO_MeanField._qtt_directed_bond_mps_from_probe_register(
        extended_result.state, register, T_U, sites, register_sites;
        cutoff=1e-12, maxdim=32,
    )
    register_up_bond = MPO_MeanField._qtt_mps_amplitudes(up_bond_register, sites)
    translations = MPO_MeanField.build_translation_square(sites)
    fock, _ = MPO_MeanField._qtt_square_fock_mpo_from_bond_mps(
        right_bond_register, up_bond_register, sites, translations, 0.7;
        cutoff=1e-12, maxdim=32,
    )
    bra_states, ket_states = MPO_MeanField.precompute_qtt_states(sites)
    for index in 0:((1 << levels) - 1)
        x, y = MPO_MeanField.square_lattice_decoder(index, levels)
        site = index + 1
        if x < (1 << (levels ÷ 2)) - 1
            right = MPO_MeanField.square_lattice_index(x + 1, y, levels)
            @test isapprox(MPO_MeanField.MatrixChecker(fock, sites, site, right, bra_states, ket_states),
                -0.7 * register_right_bond[index + 1]; atol=3e-10, rtol=3e-10)
        end
        if y < (1 << (levels ÷ 2)) - 1
            up = MPO_MeanField.square_lattice_index(x, y + 1, levels)
            @test isapprox(MPO_MeanField.MatrixChecker(fock, sites, site, up, bra_states, ket_states),
                -0.7 * register_up_bond[index + 1]; atol=3e-10, rtol=3e-10)
        end
    end

    # Register moments are the same stochastic-trace moments as an explicit
    # Walsh matrix, but retain only two QTT recurrence states.
    moment_degree = 8
    register_moments = MPO_MeanField._qtt_probe_register_moments(
        extended_H, register, moment_degree; cutoff=1e-12, maxdim=32,
    )
    exact_moments = [sum(cos(order * acos(value)) for value in diagonal)
        for order in 0:moment_degree]
    @test isapprox(register_moments, exact_moments; atol=3e-10, rtol=3e-10)
end
