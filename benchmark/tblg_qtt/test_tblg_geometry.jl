using Test
using LinearAlgebra

include(joinpath(@__DIR__, "tblg_geometry.jl"))
using .TBLGCircularFlake

@testset "TBLG circular-flake geometry" begin
    p = TBLGParameters(angle_deg=5.0)
    sites = circular_flake(20.0, p)
    @test !isempty(sites)
    @test length(unique((s.layer, s.sublattice, s.n1, s.n2) for s in sites)) == length(sites)
    @test all(hypot(s.x, s.y) <= 20.0 + 1e-10 for s in sites)
    @test opnorm(hopping_matrix(sites, p) - hopping_matrix(sites, p)', Inf) < 1e-12
    @test hopping_matrix(sites, p) ==
          TBLGCircularFlake._hopping_matrix_bruteforce(sites, p)

    # An offset centre removes the centered-shell degeneracy and reaches an
    # exact 2^10-site address space for the 5-degree test.
    radius = find_exact_radius(1024, p; radius_max=40.0, center=(0.1, 0.0))
    @test radius !== nothing
    exact_sites = circular_flake(radius, p; center=(0.1, 0.0))
    @test length(exact_sites) == 2^10
    @test validate_site_addressing(exact_sites).qtt_bits == 10

    # Morton ordering changes only the address order, not the physical site set.
    morton_sites = morton_ordered_sites(exact_sites)
    @test validate_site_addressing(morton_sites).qtt_bits == 10
    @test Set((s.layer, s.sublattice, s.n1, s.n2) for s in morton_sites) ==
          Set((s.layer, s.sublattice, s.n1, s.n2) for s in exact_sites)
    @test [s.id for s in morton_sites] == collect(1:2^10)

    # Aligned-layer vertical hopping is exactly V_pp_sigma at r=d0.
    p_aligned = TBLGParameters(angle_deg=0.0)
    aligned = circular_flake(3.0, p_aligned)
    a = findfirst(s -> s.layer == 1 && s.sublattice == 1 && s.n1 == 0 && s.n2 == 0, aligned)
    b = findfirst(s -> s.layer == 2 && s.sublattice == 1 && s.n1 == 0 && s.n2 == 0, aligned)
    @test a !== nothing && b !== nothing
    @test hopping_value(aligned[a], aligned[b], p_aligned) ≈ p_aligned.vpp_sigma_0
end

println("TBLG geometry tests passed.")
