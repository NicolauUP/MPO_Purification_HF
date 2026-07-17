function square_site_densities(sys)
    N = 2^sys.params.L
    return [
        real(MatrixChecker(sys.ρ, sys.sites, site, site, sys.bra_states, sys.ket_states))
        for site in 1:N
    ]
end

"""
    square_staggered_density(density, L)

Return the checkerboard density order parameter
``sum_{x,y} (-1)^(x+y) (n_xy - n_bar) / N`` using the package's interleaved
square-site convention. A uniform density has zero order parameter.
"""
function square_staggered_density(density, L)
    N = length(density)
    N == 2^L || throw(ArgumentError("density length must equal 2^L"))
    nbar = sum(density) / N
    return sum(begin
        x, y = square_lattice_decoder(site - 1, L)
        (-1)^(x + y) * (density[site] - nbar)
    end for site in 1:N) / N
end

@testset "P2.5 square SCF checkerboard-seed response" begin
    function checkerboard_seed(sign)
        return (x, y) -> sign * (iseven(Int(x) + Int(y)) ? 0.20 : -0.20)
    end

    # W is deliberately absent: any staggered density is selected by the
    # initial HF field, not imposed by an external staggered potential.
    function seeded_square_system(sign)
        return System(parameters_square(
            L=4,
            t=(-1.0, -1.0),
            W=nothing,
            U=3.0,
            S=checkerboard_seed(sign),
            density=0.5,
            purification_steps=45,
            itensors_maxdim=96,
            scf_mixing=0.30,
            scf_tol=0.1,
            scf_max_iterations=40,
        ))
    end

    positive = seeded_square_system(1.0)
    negative = seeded_square_system(-1.0)
    @test run_scf!(positive, -10.0, 10.0; purification_method=:sp2, verbose=:nothing)
    @test run_scf!(negative, -10.0, 10.0; purification_method=:sp2, verbose=:nothing)

    density_positive = square_site_densities(positive)
    density_negative = square_site_densities(negative)
    order_positive = square_staggered_density(density_positive, positive.params.L)
    order_negative = square_staggered_density(density_negative, negative.params.L)

    @test isapprox(sum(density_positive), 8.0; atol=3e-3)
    @test isapprox(sum(density_negative), 8.0; atol=3e-3)
    @test abs(order_positive) > 0.05
    @test abs(order_negative) > 0.05
    @test order_positive * order_negative < 0.0
end
