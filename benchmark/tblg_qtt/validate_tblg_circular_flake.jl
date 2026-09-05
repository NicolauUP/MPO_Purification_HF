#!/usr/bin/env julia

using LinearAlgebra
using Printf
using SparseArrays

include(joinpath(@__DIR__, "tblg_geometry.jl"))
using .TBLGCircularFlake

function usage()
    println("usage: julia validate_tblg_circular_flake.jl [angle_deg] [radius] [center_x] [center_y]")
    println("default: angle_deg=5.0 radius=20.0 center=(0.0,0.0)")
end

length(ARGS) > 4 && (usage(); exit(2))
angle = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 5.0
radius = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 20.0
center = (length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0,
          length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.0)
p = TBLGParameters(angle_deg=angle)
sites = circular_flake(radius, p; center=center)
counts = count_by_layer_sublattice(sites)
println("TBLG circular-flake geometry validation")
@printf("angle_deg=%.8f radius=%.8f center=(%.8f,%.8f) a0=%.8f d0=%.8f\n",
        angle, radius, center[1], center[2], p.a0, p.d0)
@printf("physical_sites=%d power_of_two=%s\n", length(sites), ispow2(length(sites)))
println("counts_by_layer_sublattice=", counts)
if ispow2(length(sites))
    println("addressing=", validate_site_addressing(sites))
else
    next_n = 2^ceil(Int, log2(length(sites)))
    @printf("addressing=not_exact_power_of_two next_quantics_dimension=%d padding_slots=%d\n",
            next_n, next_n - length(sites))
end

H = hopping_matrix(sites, p)
hermiticity = opnorm(H - H', Inf)
@printf("matrix_dimension=%d nnz=%d hermiticity_inf_norm=%.3e\n",
        size(H, 1), nnz(H), hermiticity)
hermiticity <= 1e-12 || error("hopping matrix is not Hermitian")
all(iszero, diag(H)) || error("hopping diagonal must be zero in this prototype")
