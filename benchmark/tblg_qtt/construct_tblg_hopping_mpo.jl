#!/usr/bin/env julia

"""Construct and validate a QTT MPO for a small atomistic TBLG flake.

This is an *exact-tensor* prototype: the sparse one-particle hopping matrix is
first materialized as a 2L-index ITensor and then factorized into an MPO. It
therefore tests the representation and the importance of site ordering; it is
not yet a scalable TCI construction.

Usage:
  julia --startup-file=no --project=. construct_tblg_hopping_mpo.jl \
      [target=1024] [angle_deg=5.0] [cutoff=1e-12] [maxdim=1024] [output_dir]
"""

using ITensors
using ITensorMPS
using LinearAlgebra
using Printf
using Random
using SparseArrays
using Statistics

include(joinpath(@__DIR__, "tblg_geometry.jl"))
using .TBLGCircularFlake

target = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1024
angle = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 5.0
cutoff = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-12
maxdim = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : target
output = length(ARGS) >= 5 ? abspath(ARGS[5]) :
         joinpath(@__DIR__, "results", "tblg_hopping_mpo_$(target)")

ispow2(target) || error("target must be an exact power of two, got $target")
target <= 2048 || error("this dense-tensor prototype is deliberately restricted to target <= 2048")
cutoff >= 0 || error("cutoff must be nonnegative")
maxdim > 0 || error("maxdim must be positive")
isdir(output) && error("refusing to overwrite existing output directory: $output")
mkpath(output)

"""Return a reverse Cuthill--McKee permutation for the undirected sparsity graph."""
function reverse_cuthill_mckee(H::SparseMatrixCSC)
    n = size(H, 1)
    neighbours = Vector{Vector{Int}}(undef, n)
    degrees = zeros(Int, n)
    for vertex in 1:n
        ns = findnz(H[:, vertex])[1]
        filter!(neighbour -> neighbour != vertex, ns)
        neighbours[vertex] = ns
        degrees[vertex] = length(ns)
    end

    visited = falses(n)
    order = Int[]
    while length(order) < n
        roots = findall(!, visited)
        root = roots[argmin(degrees[roots])]
        queue = [root]
        visited[root] = true
        head = 1
        while head <= length(queue)
            vertex = queue[head]
            head += 1
            push!(order, vertex)
            candidates = [neighbour for neighbour in neighbours[vertex] if !visited[neighbour]]
            sort!(candidates, by=neighbour -> (degrees[neighbour], neighbour))
            for neighbour in candidates
                visited[neighbour] && continue
                visited[neighbour] = true
                push!(queue, neighbour)
            end
        end
    end
    return reverse(order)
end

function reordered_sites(sites, permutation)
    return [SiteRecord(address, s.layer, s.sublattice, s.n1, s.n2, s.x, s.y, s.z)
            for (address, s) in enumerate(sites[permutation])]
end

function write_address_table(path, sites)
    open(path, "w") do io
        println(io, "binary_address,layer,sublattice,n1,n2,x,y,z")
        for site in sites
            println(io, join((site.id, site.layer, site.sublattice, site.n1,
                              site.n2, site.x, site.y, site.z), ','))
        end
    end
end

function binary_matrix_tensor(H::SparseMatrixCSC, sites::Vector{<:Index})
    n = size(H, 1)
    n == 1 << length(sites) || error("matrix dimension and QTT bits disagree")
    # The first L indices are bra/output bits; the last L are ket/input bits.
    tensor = ITensor(prime.(sites)..., dag.(sites)...)
    rows, cols, values = findnz(H)
    L = length(sites)
    for (row, col, value) in zip(rows, cols, values)
        bra = [(((row - 1) >> (L - bit)) & 1) + 1 for bit in 1:L]
        ket = [(((col - 1) >> (L - bit)) & 1) + 1 for bit in 1:L]
        tensor[(bra..., ket...)...] = value
    end
    return tensor
end

function mpo_matrix_element(mpo::MPO, sites::Vector{<:Index}, row::Int, col::Int)
    L = length(sites)
    value = ITensor(1.0)
    for bit in 1:L
        shift = L - bit
        bra_bit = ((row - 1) >> shift) & 1
        ket_bit = ((col - 1) >> shift) & 1
        # This is the same MSB-first convention as MatrixChecker.
        value *= mpo[bit] * dag(state(sites[bit]', bra_bit == 0 ? "0" : "1")) *
                 state(sites[bit], ket_bit == 0 ? "0" : "1")
    end
    return scalar(value)
end

function sampled_error(mpo::MPO, H::SparseMatrixCSC, sites::Vector{<:Index}; nsamples=4096)
    rows, cols, _ = findnz(H)
    rng = MersenneTwister(510578)
    count = min(nsamples, length(rows))
    nonzero_indices = rand(rng, 1:length(rows), count)
    max_error, squared_error = 0.0, 0.0
    for index in nonzero_indices
        actual = mpo_matrix_element(mpo, sites, rows[index], cols[index])
        err = abs(actual - H[rows[index], cols[index]])
        max_error = max(max_error, err)
        squared_error += err^2
    end
    # Include deterministic random locations; nearly all are exact zeros and
    # catch spurious support introduced by an inaccurate compression.
    for _ in 1:nsamples
        row, col = rand(rng, axes(H, 1)), rand(rng, axes(H, 2))
        actual = mpo_matrix_element(mpo, sites, row, col)
        err = abs(actual - H[row, col])
        max_error = max(max_error, err)
        squared_error += err^2
    end
    return (max_abs_error=max_error,
            rms_error=sqrt(squared_error / (count + nsamples)),
            nonzero_samples=count,
            total_samples=count + nsamples)
end

function bandwidth(H)
    rows, cols, values = findnz(H)
    return isempty(values) ? 0 : maximum(abs.(rows .- cols))
end

p = TBLGParameters(angle_deg=angle)
center = (0.1, 0.0)
radius = find_exact_radius(target, p; radius_max=40.0, center=center)
radius === nothing && error("could not find a radius with exactly $target sites")
coordinate_sites = circular_flake(radius, p; center=center)
H_coordinate = hopping_matrix(coordinate_sites, p)
rcm_permutation = reverse_cuthill_mckee(H_coordinate)
rcm_sites = reordered_sites(coordinate_sites, rcm_permutation)
H_rcm = H_coordinate[rcm_permutation, rcm_permutation]
write_address_table(joinpath(output, "coordinate_binary_address.csv"), coordinate_sites)
write_address_table(joinpath(output, "rcm_binary_address.csv"), rcm_sites)

sites = siteinds("Qubit", quantics_level(target))
# The dense prototype starts from one 2L-index tensor. Its high intermediate
# order is intentional and bounded by target <= 2048, so suppress diagnostic
# warnings that otherwise obscure the benchmark output.
ITensors.set_warn_order(2length(sites) + 1)
summary_rows = String[]
for (name, H) in (("coordinate", H_coordinate), ("rcm", H_rcm))
    record = nothing
    elapsed = @elapsed begin
        tensor = binary_matrix_tensor(H, sites)
        mpo = MPO(tensor, sites; cutoff=cutoff, maxdim=maxdim)
        validation = sampled_error(mpo, H, sites)
        dimensions = linkdims(mpo)
        record = (max_chi=maxlinkdim(mpo), mean_chi=mean(dimensions),
                  validation=validation, dimensions=dimensions)
        open(joinpath(output, "$(name)_bond_dimensions.csv"), "w") do bonds
            println(bonds, "bond,dimension")
            for (bond, dimension) in enumerate(dimensions)
                println(bonds, "$(bond),$(dimension)")
            end
        end
    end
    push!(summary_rows, @sprintf("%s,%d,%d,%.12g,%.12e,%.12e,%d,%d,%.9f",
        name, bandwidth(H), record.max_chi, record.mean_chi,
        record.validation.max_abs_error, record.validation.rms_error,
        record.validation.nonzero_samples, record.validation.total_samples, elapsed))
end
open(joinpath(output, "summary.csv"), "w") do io
    println(io, "ordering,bandwidth,max_chi,mean_chi,sampled_max_abs_error,sampled_rms_error,nonzero_samples,total_samples,construction_time_s")
    foreach(row -> println(io, row), summary_rows)
end

open(joinpath(output, "metadata.toml"), "w") do io
    println(io, "target_sites = $target")
    println(io, "qtt_bits = $(quantics_level(target))")
    println(io, "angle_deg = $angle")
    println(io, "radius = $radius")
    println(io, "center = [$(center[1]), $(center[2])]")
    println(io, "cutoff = $cutoff")
    println(io, "maxdim = $maxdim")
    println(io, "nnz = $(nnz(H_coordinate))")
    println(io, "hermiticity_residual = $(opnorm(H_coordinate - H_coordinate', Inf))")
    println(io, "construction = \"dense 2L-index tensor followed by MPO SVD; prototype only\"")
    println(io, "validation = \"4096 sampled nonzero entries plus 4096 deterministic random entries\"")
end

println("output=$output")
println(read(joinpath(output, "summary.csv"), String))
