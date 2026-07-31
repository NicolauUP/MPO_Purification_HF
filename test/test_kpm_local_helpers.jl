using Test
using LinearAlgebra
using Random
using Statistics

@testset "Jackson KPM helper conventions" begin
    include(joinpath(
        @__DIR__, "..", "runs", "common", "kpm_local_helpers.jl",
    ))

    function allocating_kpm_apply(H, probes, coefficients)
        previous = copy(probes)
        current = H * probes
        result = (coefficients[1] / 2) .* previous .+
                 coefficients[2] .* current
        for order in 2:(length(coefficients) - 1)
            following = 2 .* (H * current) .- previous
            result .+= coefficients[order + 1] .* following
            previous, current = current, following
        end
        return result
    end

    function allocating_trace_moments(H, probes, degree)
        normalization = inv(size(probes, 2))
        moments = Vector{Float64}(undef, degree + 1)
        previous = copy(probes)
        moments[1] = real(dot(vec(probes), vec(previous))) * normalization
        degree == 0 && return moments
        current = H * probes
        moments[2] =
            real(dot(vec(probes), vec(current))) * normalization
        for order in 2:degree
            following = 2 .* (H * current) .- previous
            moments[order + 1] =
                real(dot(vec(probes), vec(following))) * normalization
            previous, current = current, following
        end
        return moments
    end

    degree = 64
    coefficients = projector_coefficients(degree, 0.0)
    polynomial(x) = begin
        previous = 1.0
        current = x
        value = coefficients[1] / 2 + coefficients[2] * current
        for order in 2:degree
            following = 2x * current - previous
            value += coefficients[order + 1] * following
            previous, current = current, following
        end
        value
    end
    @test polynomial(-0.8) > 0.99
    @test polynomial(0.8) < 0.01
    @test isapprox(polynomial(-0.4), 1 - polynomial(0.4); atol=1e-12)

    probes = probing_matrix(8, 8, :hadamard, 42)
    @test isapprox(probes * probes' / 8, Matrix{Float64}(I, 8, 8);
        atol=1e-14)
    @test probing_matrix(8, 3, :hadamard, 42) == probes[:, 1:3]

    codes = [0, 2, 4, 6, 1, 3, 5, 7]
    coded = coded_hadamard_matrix(codes, 8, 42)
    @test isapprox(coded * coded' / 8, Matrix{Float64}(I, 8, 8);
        atol=1e-14)
    @test coded_hadamard_matrix(codes, 3, 42) == coded[:, 1:3]
    @test hcat(
        coded_hadamard_block(codes, 1, 3, 42),
        coded_hadamard_block(codes, 4, 5, 42),
    ) == coded
    @test hcat(
        hadamard_block(8, 1, 3, 42),
        hadamard_block(8, 4, 5, 42),
    ) == probes

    rectangular_codes = [
        coordinate_interleaved_code(x, y, 3, 2)
        for y in 0:3 for x in 0:7
    ]
    @test sort(rectangular_codes) == collect(0:31)
    @test coordinate_interleaved_code(5, 2, 3, 2) == 22
    # The production 1024x512 ordering used before the rectangular
    # generalization interleaved y/x bits 0:8, then appended x bit 9.
    legacy_1024x512_code(x, y) = begin
        code = 0
        for bit in 0:8
            code |= ((y >> bit) & 1) << (2bit)
            code |= ((x >> bit) & 1) << (2bit + 1)
        end
        code |= ((x >> 9) & 1) << 18
        code
    end
    for x in (0, 1, 17, 1023), y in (0, 1, 19, 511)
        @test coordinate_interleaved_code(x, y, 10, 9) ==
              legacy_1024x512_code(x, y)
    end
    square_codes = [
        coordinate_interleaved_code(x, y, 2, 2)
        for y in 0:3 for x in 0:3
    ]
    @test sort(square_codes) == collect(0:15)

    block_hamiltonian = SymTridiagonal(
        collect(range(-0.7, 0.7; length=8)),
        fill(0.05, 7),
    )
    block_coefficients = projector_coefficients(32, 0.0)
    full_moments =
        kpm_trace_moments(block_hamiltonian, coded, 32)
    blocked_moments = zeros(33)
    blocked_density_numerator = zeros(8)
    for block in (
        coded_hadamard_block(codes, 1, 3, 42),
        coded_hadamard_block(codes, 4, 5, 42),
    )
        blocked_moments .+=
            size(block, 2) .* kpm_trace_moments(
                block_hamiltonian, block, 32,
            )
        blocked_density_numerator .+= vec(sum(
            kpm_apply(block_hamiltonian, block, block_coefficients) .* block;
            dims=2,
        ))
    end
    blocked_moments ./= size(coded, 2)
    full_filtered =
        kpm_apply(block_hamiltonian, coded, block_coefficients)
    full_density =
        vec(mean(full_filtered .* coded; dims=2))
    @test isapprox(blocked_moments, full_moments; atol=1e-12)
    @test isapprox(
        blocked_density_numerator ./ size(coded, 2),
        full_density;
        atol=1e-12,
    )
    @test isapprox(
        kpm_apply(block_hamiltonian, coded, block_coefficients),
        allocating_kpm_apply(
            block_hamiltonian, coded, block_coefficients,
        );
        atol=1e-13,
        rtol=1e-13,
    )
    @test isapprox(
        kpm_trace_moments(block_hamiltonian, coded, 32),
        allocating_trace_moments(block_hamiltonian, coded, 32);
        atol=1e-13,
        rtol=1e-13,
    )

    diagonal = Diagonal([-0.8, -0.3, 0.3, 0.8])
    trace_probes = probing_matrix(4, 4, :hadamard, 17)
    trace_moments = kpm_trace_moments(diagonal, trace_probes, 128)
    chemical_potential = find_scaled_chemical_potential(
        trace_moments, 2.0; tolerance=1e-10,
    )
    @test abs(chemical_potential.scaled_mu) < 1e-10
    @test isapprox(chemical_potential.trace, 2.0; atol=1e-10)
    @test isapprox(trace_moments[1], 4.0; atol=1e-14)
    @test isapprox(trace_moments[2], 0.0; atol=1e-14)
end
