@testset "QTT charge analysis" begin
    L = 4
    params = ParametersSquare(
        L=L,
        t=(-1.0, -1.0),
        U=0.0,
        W=nothing,
        S=nothing,
        tci_tol=1e-14,
        itensors_tol=1e-12,
        itensors_maxdim=64,
        density=0.5,
        purification_steps=20,
        scf_mixing=0.5,
        scf_tol=1e-8,
        scf_max_iterations=10,
    )
    system = System(params)
    density_function(z) = begin
        x, y = square_lattice_decoder(Int(z), L)
        0.4 + 0.03x + 0.02y + 0.01x * y
    end
    _, density_mpo, _ = MPO_MeanField.Quantics_TCI(
        density_function, Float64, system.sites, 1e-14,
    )
    system.ρ = density_mpo

    density = density_diagonal_mps(system)
    values = [density_function(index) for index in 0:15]
    recovered = [qtt_mps_amplitude(density, system.sites, index) for index in 0:15]
    @test recovered ≈ values atol=1e-12

    centered, mean_density, trace_value = centered_density_mps(
        density, system.sites; cutoff=1e-14, maxdim=64,
    )
    delta = values .- sum(values) / length(values)
    @test real(trace_value) ≈ sum(values) atol=1e-12
    @test real(mean_density) ≈ sum(values) / length(values) atol=1e-12

    participation = qtt_charge_ipr(centered; cutoff=1e-14, maxdim=64)
    expected_ipr = sum(abs2.(delta) .^ 2) / sum(abs2, delta)^2
    @test participation.ipr ≈ expected_ipr rtol=1e-10
    @test participation.participation ≈ inv(expected_ipr) rtol=1e-10

    real_d2 = qtt_multiscale_d2(
        centered; cutoff=1e-14, maxdim=64, keep=:prefix, fit_scales=1:2,
    )
    normalized_mass = abs2.(delta) ./ sum(abs2, delta)
    direct_real_z2 = Float64[]
    for scale in 0:2
        box_side = 1 << (2 - scale)
        masses = [sum(
            normalized_mass[square_lattice_index(x, y, L)]
            for x in bx:(bx + box_side - 1), y in by:(by + box_side - 1)
        ) for bx in 0:box_side:3, by in 0:box_side:3]
        push!(direct_real_z2, sum(abs2, masses))
    end
    @test [row.z2 for row in real_d2.scales] ≈ direct_real_z2 rtol=1e-10

    transformed = qtt_fourier_square(
        centered, system.sites;
        sign=-1, cutoff_MPO=1e-14, cutoff=1e-14, maxdim=64,
    )
    side = 4
    for kx in 0:(side - 1), ky in 0:(side - 1)
        expected = sum(
            delta[square_lattice_index(x, y, L)] *
            exp(-2im * pi * (kx * x + ky * y) / side)
            for x in 0:(side - 1), y in 0:(side - 1)
        ) / side
        @test qtt_fourier_amplitude(transformed, kx, ky) ≈ expected atol=1e-10
    end
    @test real(inner(transformed, transformed)) ≈ sum(abs2, delta) rtol=1e-10


    fourier_d2 = qtt_multiscale_d2(
        transformed; cutoff=1e-14, maxdim=64, keep=:suffix, fit_scales=1:2,
    )
    direct_fourier = [
        sum(
            delta[square_lattice_index(x, y, L)] *
            exp(-2im * pi * (kx * x + ky * y) / side)
            for x in 0:(side - 1), y in 0:(side - 1)
        ) / side
        for kx in 0:(side - 1), ky in 0:(side - 1)
    ]
    fourier_mass = abs2.(direct_fourier) ./ sum(abs2, direct_fourier)
    direct_fourier_z2 = Float64[]
    for scale in 0:2
        box_side = 1 << (2 - scale)
        masses = [sum(
            fourier_mass[kx + 1, ky + 1]
            for kx in bx:(bx + box_side - 1), ky in by:(by + box_side - 1)
        ) for bx in 0:box_side:3, by in 0:box_side:3]
        push!(direct_fourier_z2, sum(abs2, masses))
    end
    @test [row.z2 for row in fourier_d2.scales] ≈ direct_fourier_z2 rtol=1e-10
end
