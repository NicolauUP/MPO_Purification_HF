#!/usr/bin/env julia

"""Validate a spatial binary address for an exact-power-of-two TBLG flake."""

using LinearAlgebra
using Printf
using SparseArrays

include(joinpath(@__DIR__, "tblg_geometry.jl"))
using .TBLGCircularFlake

target = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1024
angle = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 5.0
center = (length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.1,
          length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.0)
p = TBLGParameters(angle_deg=angle)
radius = find_exact_radius(target, p; radius_max=40.0, center=center)
radius === nothing && error("no exact radius found for target=$target")
sites = circular_flake(radius, p; center=center)
morton = morton_ordered_sites(sites)
validate_site_addressing(morton)

key(s) = (s.layer, s.sublattice, s.n1, s.n2)
original_by_key = Dict(key(s) => s.id for s in sites)
morton_by_key = Dict(key(s) => s.id for s in morton)
Set(keys(original_by_key)) == Set(keys(morton_by_key)) || error("Morton order changed site set")

H = hopping_matrix(sites, p)
Hm = hopping_matrix(morton, p)
findnz_bandwidth(matrix) = begin
    rows, cols, values = findnz(matrix)
    isempty(values) ? 0 : maximum(abs.(rows .- cols))
end
original_bandwidth = findnz_bandwidth(H)
morton_bandwidth = findnz_bandwidth(Hm)

function permutation_error(H, Hm, sites, morton_by_key, key)
    max_error = 0.0
    # Check that the Morton matrix is exactly the same operator under
    # permutation of the physical site addresses.
    for (i, j, value) in zip(findnz(H)...)
        mi, mj = morton_by_key[key(sites[i])], morton_by_key[key(sites[j])]
        max_error = max(max_error, abs(value - Hm[mi, mj]))
    end
    return max_error
end
max_permuted_error = permutation_error(H, Hm, sites, morton_by_key, key)

@printf("target=%d angle_deg=%.6f radius=%.15f center=(%.6f,%.6f)\n",
        target, angle, radius, center[1], center[2])
@printf("sites=%d qtt_bits=%d nnz=%d hermiticity=%.3e\n",
        length(morton), quantics_level(length(morton)), nnz(Hm), opnorm(Hm - Hm', Inf))
@printf("bandwidth_spatial=%d bandwidth_morton=%d reduction=%.3f\n",
        original_bandwidth, morton_bandwidth,
        original_bandwidth == 0 ? 0.0 : 1 - morton_bandwidth / original_bandwidth)
@printf("permuted_operator_max_error=%.3e\n", max_permuted_error)
max_permuted_error <= 1e-12 || error("Morton permutation changed hopping operator")

output = length(ARGS) >= 5 ? abspath(ARGS[5]) : joinpath(@__DIR__, "tblg_binary_address_$(target).csv")
open(output, "w") do io
    println(io, "binary_address,layer,sublattice,n1,n2,x,y,z")
    for s in morton
        println(io, join((s.id, s.layer, s.sublattice, s.n1, s.n2, s.x, s.y, s.z), ','))
    end
end
println("address_table=$output")
