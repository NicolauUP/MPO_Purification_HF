#!/usr/bin/env julia

"""Cross-interpolate an atomistic TBLG hopping oracle directly as a QTT MPO.

The row and column binary digits are fused at each QTT level. TCI therefore
constructs a rank-(left, row-bit, column-bit, right) tensor train, which maps
directly to an MPO tensor with a primed output and unprimed input Qubit index.
Unlike construct_tblg_hopping_mpo.jl, this path never materializes the dense
N by N hopping matrix during construction.
"""

using ITensors
using ITensorMPS
using LinearAlgebra
using Printf
using QuanticsTCI
using Random
using SparseArrays
import TensorCrossInterpolation as TCI

include(joinpath(@__DIR__, "tblg_geometry.jl"))
using .TBLGCircularFlake

target = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1024
angle = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 5.0
tolerance = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 1e-10
maxdim = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 512
output = length(ARGS) >= 5 ? abspath(ARGS[5]) :
         joinpath(@__DIR__, "results", "tblg_hopping_tci_$(target)")
seed = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 510578
pivot_mode = length(ARGS) >= 7 ? Symbol(ARGS[7]) : :random
radius_max = length(ARGS) >= 8 ? parse(Float64, ARGS[8]) : max(40.0, 0.8 * sqrt(target))
component = length(ARGS) >= 9 ? Symbol(ARGS[9]) : :full
ordering = length(ARGS) >= 10 ? Symbol(ARGS[10]) : :coordinate

ispow2(target) || error("target must be an exact power of two, got $target")
tolerance > 0 || error("tolerance must be positive")
maxdim > 0 || error("maxdim must be positive")
pivot_mode in (:random, :bond_multiscale) ||
    error("pivot_mode must be random or bond_multiscale")
radius_max > 0 || error("radius_max must be positive")
component in (:full, :interlayer) || error("component must be full or interlayer")
ordering in (:coordinate, :rcm) || error("ordering must be coordinate or rcm")
isdir(output) && error("refusing to overwrite existing output directory: $output")
mkpath(output)

function mpo_matrix_element(mpo::MPO, sites::Vector{<:Index}, row::Int, col::Int)
    L = length(sites)
    value = ITensor(1.0)
    for bit in 1:L
        shift = L - bit
        bra_bit = ((row - 1) >> shift) & 1
        ket_bit = ((col - 1) >> shift) & 1
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
    maximum_error, squared_error = 0.0, 0.0
    for index in nonzero_indices
        err = abs(mpo_matrix_element(mpo, sites, rows[index], cols[index]) - H[rows[index], cols[index]])
        maximum_error = max(maximum_error, err)
        squared_error += err^2
    end
    for _ in 1:nsamples
        row, col = rand(rng, axes(H, 1)), rand(rng, axes(H, 2))
        err = abs(mpo_matrix_element(mpo, sites, row, col) - H[row, col])
        maximum_error = max(maximum_error, err)
        squared_error += err^2
    end
    return (max_abs_error=maximum_error,
            rms_error=sqrt(squared_error / (count + nsamples)),
            nonzero_samples=count,
            total_samples=count + nsamples)
end

function sampled_qtt_error(qtt, H::SparseMatrixCSC; nsamples=4096)
    rows, cols, _ = findnz(H)
    rng = MersenneTwister(510578)
    count = min(nsamples, length(rows))
    nonzero_indices = rand(rng, 1:length(rows), count)
    maximum_error, squared_error = 0.0, 0.0
    for index in nonzero_indices
        err = abs(qtt(rows[index], cols[index]) - H[rows[index], cols[index]])
        maximum_error = max(maximum_error, err)
        squared_error += err^2
    end
    for _ in 1:nsamples
        row, col = rand(rng, axes(H, 1)), rand(rng, axes(H, 2))
        err = abs(qtt(row, col) - H[row, col])
        maximum_error = max(maximum_error, err)
        squared_error += err^2
    end
    return (max_abs_error=maximum_error,
            rms_error=sqrt(squared_error / (count + nsamples)))
end

"""
    bond_multiscale_pivots(H, atom_sites, L, seed; per_scale_budget=32)

Select physical nonzero hopping pairs for TCI initialization. For every binary
prefix length, representatives are drawn separately from intralayer and
interlayer bonds. Both orientations are retained, so the initial skeleton does
not break Hermiticity by construction. This sparse graph is used only for
pivot selection and the later audit, never as the tensor TCI factorizes.
"""
function bond_multiscale_pivots(
    H::SparseMatrixCSC,
    atom_sites,
    L::Int,
    seed::Int;
    per_scale_budget::Int=32,
)
    rows, cols, _ = findnz(H)
    rng = MersenneTwister(seed)
    selected = Tuple{Int,Int}[]
    seen = Set{Tuple{Int,Int}}()
    function add_pair!(row, col)
        for pair in ((row, col), (col, row))
            pair in seen && continue
            push!(seen, pair)
            push!(selected, pair)
        end
    end

    for level in 1:L
        shift = L - level
        for same_layer in (true, false)
            representatives = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
            for index in eachindex(rows)
                row, col = rows[index], cols[index]
                (atom_sites[row].layer == atom_sites[col].layer) == same_layer || continue
                key = ((row - 1) >> shift, (col - 1) >> shift)
                haskey(representatives, key) ||
                    (representatives[key] = (row, col))
            end
            keys_at_scale = collect(keys(representatives))
            shuffle!(rng, keys_at_scale)
            for key in Iterators.take(keys_at_scale, per_scale_budget)
                add_pair!(representatives[key]...)
            end
        end
    end
    return [Int[row, col] for (row, col) in selected]
end

p = TBLGParameters(angle_deg=angle)
center = (0.1, 0.0)
radius = find_exact_radius(target, p; radius_max=radius_max, center=center)
radius === nothing && error("could not find a radius with exactly $target sites")
atom_sites = circular_flake(radius, p; center=center)
validate_site_addressing(atom_sites)

# The oracle is the only input to cross interpolation. QuanticsTCI uses
# one-based grid coordinates, precisely matching the site addresses.
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
            vertex = queue[head]; head += 1; push!(order, vertex)
            candidates = [x for x in neighbours[vertex] if !visited[x]]
            sort!(candidates, by=x -> (degrees[x], x))
            for x in candidates
                visited[x] && continue
                visited[x] = true
                push!(queue, x)
            end
        end
    end
    return reverse(order)
end
reference = pivot_mode === :bond_multiscale ?
    hopping_matrix(atom_sites, p; include_intralayer=component === :full,
                   include_interlayer=true) : nothing
if ordering === :rcm
    reference === nothing && (reference = hopping_matrix(
        atom_sites, p; include_intralayer=component === :full,
        include_interlayer=true))
    permutation = reverse_cuthill_mckee(reference)
    atom_sites = [SiteRecord(address, s.layer, s.sublattice, s.n1, s.n2,
                             s.x, s.y, s.z)
                  for (address, s) in enumerate(atom_sites[permutation])]
    reference = reference[permutation, permutation]
end
hopping_oracle(row::Int, col::Int) = begin
    a, b = atom_sites[row], atom_sites[col]
    component === :interlayer && a.layer == b.layer ? 0.0 : hopping_value(a, b, p)
end
initial_pivots = pivot_mode === :bond_multiscale ?
    bond_multiscale_pivots(reference, atom_sites, quantics_level(target), seed) : nothing
qtt = nothing
ranks = nothing
errors = nothing
Random.seed!(seed)
elapsed = @elapsed begin
    qtt, ranks, errors = quanticscrossinterpolate(
        Float64, hopping_oracle, (target, target), initial_pivots;
        unfoldingscheme=:fused,
        tolerance=tolerance,
        maxbonddim=maxdim,
    )
end

tt_fused = TCI.tensortrain(qtt.tci)
L = quantics_level(target)
# Fused row/column local dimension 4 is reshaped with row first, because the
# Quantics fused grid makes its first coordinate vary fastest.
tt_mpo = TCI.TensorTrain{4}(tt_fused, fill([2, 2], L))
qtt_sites = siteinds("Qubit", L)
hopping_mpo = MPO(tt_mpo; sites=[[prime(site), dag(site)] for site in qtt_sites])

# This sparse matrix is constructed only after TCI, solely for the small-system
# audit. The interpolation itself never sees it.
reference === nothing && (reference = hopping_matrix(
    atom_sites, p;
    include_intralayer=component === :full,
    include_interlayer=true,
))
validation = sampled_error(hopping_mpo, reference, qtt_sites)
qtt_validation = sampled_qtt_error(qtt, reference)

open(joinpath(output, "summary.toml"), "w") do io
    println(io, "target_sites = $target")
    println(io, "qtt_bits = $L")
    println(io, "angle_deg = $angle")
    println(io, "radius = $radius")
    println(io, "radius_search_max = $radius_max")
    println(io, "component = \"$component\"")
    println(io, "ordering = \"$ordering\"")
    println(io, "tci_tolerance = $tolerance")
    println(io, "requested_maxdim = $maxdim")
    println(io, "seed = $seed")
    println(io, "pivot_mode = \"$pivot_mode\"")
    println(io, "initial_pivot_count = $(isnothing(initial_pivots) ? 0 : length(initial_pivots))")
    println(io, "tci_max_chi = $(maximum(TCI.linkdims(tt_mpo); init=1))")
    println(io, "mpo_max_chi = $(maxlinkdim(hopping_mpo))")
    println(io, "mpo_mean_chi = $(sum(linkdims(hopping_mpo)) / length(linkdims(hopping_mpo)))")
    println(io, "tci_construction_time_s = $elapsed")
    println(io, "sampled_max_abs_error = $(validation.max_abs_error)")
    println(io, "sampled_rms_error = $(validation.rms_error)")
    println(io, "qtt_sampled_max_abs_error = $(qtt_validation.max_abs_error)")
    println(io, "qtt_sampled_rms_error = $(qtt_validation.rms_error)")
    println(io, "nonzero_samples = $(validation.nonzero_samples)")
    println(io, "total_samples = $(validation.total_samples)")
    println(io, "last_reported_tci_error = $(isempty(errors) ? NaN : errors[end])")
    println(io, "construction = \"fused row/column QuanticsTCI hopping oracle; no dense matrix during TCI\"")
end
open(joinpath(output, "bond_dimensions.csv"), "w") do io
    println(io, "bond,dimension")
    for (bond, dimension) in enumerate(linkdims(hopping_mpo))
        println(io, "$bond,$dimension")
    end
end
open(joinpath(output, "tci_history.csv"), "w") do io
    println(io, "step,rank,error")
    for step in eachindex(ranks)
        println(io, "$(step),$(ranks[step]),$(errors[step])")
    end
end

println("output=$output")
println(read(joinpath(output, "summary.toml"), String))
