@testset "P1.2 square Hartree dense reference" begin
    params = parameters_square(L=4, U=0.7, S=nothing)
    sys = System(params)
    occupations = [0.05 + 0.03 * i for i in 1:16]
    sys.ρ = MPO_MeanField.diagonal_mpo_from_function(
        x -> occupations[Int(x) + 1], Float64, sys.sites, params.tci_tol,
    )

    hartree = dense_matrix(extract_hartree_mpo_square(sys), sys)
    expected = zeros(16, 16)
    for site in 1:16
        for neighbour in values(square_neighbours(site, params.L))
            !isnothing(neighbour) && (expected[site, site] += params.U * occupations[neighbour])
        end
    end
    @test opnorm(hartree - expected) < 1e-10
    @test opnorm(hartree - hartree') < 1e-12

    evaluator = MPO_MeanField.HartreeEvaluateSquare(sys, Dict{Int,Float64}())
    evaluator(0.0)
    evaluator(1.0)
    expected_cached_sites = Int[
        neighbour for site in (1, 2)
        for neighbour in values(square_neighbours(site, params.L))
        if !isnothing(neighbour)
    ]
    @test sort(collect(keys(evaluator.density_cache))) == sort(unique(expected_cached_sites))
end
