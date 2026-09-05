module TBLGCircularFlake

using LinearAlgebra
using SparseArrays

export TBLGParameters, SiteRecord, circular_flake, hopping_value,
       hopping_matrix, quantics_level, validate_site_addressing,
       count_by_layer_sublattice, find_exact_radius, morton_ordered_sites,
       site_binary_address

"""
Parameters for the exploratory atomistic TBLG circular-flake prototype.

a0 is the carbon-carbon distance and d0 the layer separation. The
interlayer hopping uses an exponential Slater-Koster form. These defaults are
explicit prototype values, not parameters of the existing square model.
"""
Base.@kwdef struct TBLGParameters
    angle_deg::Float64 = 5.0
    a0::Float64 = 1.42
    d0::Float64 = 3.35
    intralayer_hopping::Float64 = -2.70
    vpp_pi_0::Float64 = -2.70
    vpp_sigma_0::Float64 = 0.48
    decay_length::Float64 = 0.319 * a0
    interlayer_cutoff::Float64 = 6.0
end

struct SiteRecord
    id::Int
    layer::Int
    sublattice::Int
    n1::Int
    n2::Int
    x::Float64
    y::Float64
    z::Float64
end

rotation_matrix(angle_deg) = begin
    θ = deg2rad(angle_deg)
    c, s = cos(θ), sin(θ)
    [c -s; s c]
end

function _layer_sites(layer::Int, angle_deg::Float64, radius::Float64,
                      p::TBLGParameters, center::Tuple{Float64,Float64})
    # Primitive vectors for a honeycomb lattice with nearest-neighbour distance a0.
    a1 = p.a0 .* [sqrt(3) / 2, 3 / 2]
    a2 = p.a0 .* [-sqrt(3) / 2, 3 / 2]
    basis = ([0.0, 0.0], [0.0, p.a0]) # A=1, B=2
    rotation = rotation_matrix(angle_deg)
    nmax = ceil(Int, radius / p.a0) + 4
    sites = SiteRecord[]
    for n1 in -nmax:nmax, n2 in -nmax:nmax, sublattice in 1:2
        r = n1 .* a1 .+ n2 .* a2 .+ basis[sublattice]
        r = rotation * r
        hypot(r[1] - center[1], r[2] - center[2]) <= radius + 10eps(radius) || continue
        push!(sites, SiteRecord(0, layer, sublattice, n1, n2, r[1], r[2],
                                layer == 1 ? 0.0 : p.d0))
    end
    return sites
end

"""
    circular_flake(radius, p)

Generate an unrelaxed two-layer TBLG circular flake. Layer 1 is unrotated and
layer 2 is rotated by p.angle_deg around the common origin. The returned
records are sorted by spatial coordinates, then layer, sublattice, and lattice
coordinates. The id field is assigned after sorting.
"""
function circular_flake(radius::Real, p::TBLGParameters=TBLGParameters();
                        center::Tuple{<:Real,<:Real}=(0.0, 0.0))
    radius > 0 || throw(ArgumentError("radius must be positive"))
    c = (Float64(center[1]), Float64(center[2]))
    raw = vcat(_layer_sites(1, 0.0, Float64(radius), p, c),
               _layer_sites(2, p.angle_deg, Float64(radius), p, c))
    # Deterministic spatial ordering. A true Morton/Hilbert key will be added
    # after the geometry/count validation passes.
    sort!(raw, by=s -> (round(s.x, digits=10), round(s.y, digits=10),
                        s.layer, s.sublattice, s.n1, s.n2))
    return [SiteRecord(i, s.layer, s.sublattice, s.n1, s.n2, s.x, s.y, s.z)
            for (i, s) in enumerate(raw)]
end

"""Interleave the low-order bits of two nonnegative integers (Morton key)."""
function _interleave_bits(x::Integer, y::Integer, bits::Integer)
    key = 0
    for bit in 0:(bits - 1)
        key |= ((x >> bit) & 1) << (2bit + 1)
        key |= ((y >> bit) & 1) << (2bit)
    end
    return key
end

"""
    morton_ordered_sites(sites; coordinate_bits=16)

Return the same physical sites in a deterministic spatial Morton order and
reassign one-based binary addresses. Coordinates are quantized only for
ordering; layer, sublattice, and lattice coordinates break coincident-key ties.
The resulting address is the QTT index `address-1`, and therefore requires
`log2(length(sites))` bits when the site count is a power of two.
"""
function morton_ordered_sites(sites::AbstractVector{<:SiteRecord};
                              coordinate_bits::Integer=16)
    isempty(sites) && return SiteRecord[]
    coordinate_bits > 0 || throw(ArgumentError("coordinate_bits must be positive"))
    xmin, xmax = extrema(s.x for s in sites)
    ymin, ymax = extrema(s.y for s in sites)
    levels = (1 << coordinate_bits) - 1
    quantize(value, lo, hi) = hi == lo ? 0 : clamp(round(Int, levels * (value - lo) / (hi - lo)), 0, levels)
    keyed = [(_interleave_bits(quantize(s.x, xmin, xmax), quantize(s.y, ymin, ymax), coordinate_bits),
              s.layer, s.sublattice, s.n1, s.n2, s) for s in sites]
    sort!(keyed, by=k -> k[1:5])
    return [SiteRecord(i, s.layer, s.sublattice, s.n1, s.n2, s.x, s.y, s.z)
            for (i, (_, _, _, _, _, s)) in enumerate(keyed)]
end

"""Return the one-based binary/QTT address of a site id in an ordered list."""
function site_binary_address(sites::AbstractVector{<:SiteRecord}, site_id::Integer)
    1 <= site_id <= length(sites) || throw(BoundsError(sites, site_id))
    return sites[site_id].id
end

function _same_intralayer_bond(a::SiteRecord, b::SiteRecord)
    a.layer == b.layer || return false
    if a.sublattice == 1 && b.sublattice == 2
        return (b.n1 == a.n1 && b.n2 == a.n2) ||
               (b.n1 == a.n1 - 1 && b.n2 == a.n2) ||
               (b.n1 == a.n1 && b.n2 == a.n2 - 1)
    elseif a.sublattice == 2 && b.sublattice == 1
        return _same_intralayer_bond(b, a)
    end
    return false
end

function _interlayer_hopping(a::SiteRecord, b::SiteRecord, p::TBLGParameters)
    a.layer == b.layer && return 0.0
    dx, dy, dz = a.x - b.x, a.y - b.y, a.z - b.z
    r = sqrt(dx^2 + dy^2 + dz^2)
    r <= p.interlayer_cutoff || return 0.0
    c2 = (dz / r)^2
    vpi = p.vpp_pi_0 * exp(-(r - p.a0) / p.decay_length)
    vsigma = p.vpp_sigma_0 * exp(-(r - p.d0) / p.decay_length)
    return vpi * (1 - c2) + vsigma * c2
end

"""
    hopping_value(a, b, p)

Return the real single-particle hopping between two site records. Intralayer
terms are nearest-neighbour p.intralayer_hopping; interlayer terms use the
distance/angular Slater-Koster law and cutoff.
"""
function hopping_value(a::SiteRecord, b::SiteRecord, p::TBLGParameters=TBLGParameters())
    a.id == b.id && return 0.0
    if a.layer == b.layer
        return _same_intralayer_bond(a, b) ? p.intralayer_hopping : 0.0
    end
    return _interlayer_hopping(a, b, p)
end

function hopping_matrix(sites::AbstractVector{<:SiteRecord},
                        p::TBLGParameters=TBLGParameters();
                        include_intralayer::Bool=true,
                        include_interlayer::Bool=true)
    n = length(sites)
    rows, cols, vals = Int[], Int[], Float64[]
    # Intralayer nearest-neighbour bonds are known exactly from the honeycomb
    # lattice coordinates. Iterating only A sites avoids duplicate bonds.
    addresses = Dict{Tuple{Int,Int,Int,Int},Int}()
    for site in sites
        include_intralayer || break
        addresses[(site.layer, site.sublattice, site.n1, site.n2)] = site.id
    end
    for site in sites
        site.sublattice == 1 || continue
        for (n1, n2) in ((site.n1, site.n2),
                         (site.n1 - 1, site.n2),
                         (site.n1, site.n2 - 1))
            key = (site.layer, 2, n1, n2)
            haskey(addresses, key) || continue
            neighbour = addresses[key]
            push!(rows, site.id); push!(cols, neighbour); push!(vals, p.intralayer_hopping)
            push!(rows, neighbour); push!(cols, site.id); push!(vals, p.intralayer_hopping)
        end
    end

    # The interlayer cutoff makes this a local geometric neighbour problem,
    # not an all-pairs problem. Bin layer 2 in cells of the cutoff width and
    # inspect only the nine neighbouring cells for each layer-1 site.
    include_interlayer || return sparse(rows, cols, vals, n, n)
    cell_width = p.interlayer_cutoff
    layer2_cells = Dict{Tuple{Int,Int},Vector{Int}}()
    for site in sites
        site.layer == 2 || continue
        key = (floor(Int, site.x / cell_width), floor(Int, site.y / cell_width))
        push!(get!(layer2_cells, key, Int[]), site.id)
    end
    for site in sites
        site.layer == 1 || continue
        cx, cy = floor(Int, site.x / cell_width), floor(Int, site.y / cell_width)
        for dx in -1:1, dy in -1:1
            for neighbour in get(layer2_cells, (cx + dx, cy + dy), Int[])
                t = _interlayer_hopping(site, sites[neighbour], p)
                iszero(t) && continue
                push!(rows, site.id); push!(cols, neighbour); push!(vals, t)
                push!(rows, neighbour); push!(cols, site.id); push!(vals, t)
            end
        end
    end
    return sparse(rows, cols, vals, n, n)
end

# Retained solely for small-system regression checks of the local neighbour
# construction above. It is intentionally not part of the public API.
function _hopping_matrix_bruteforce(sites::AbstractVector{<:SiteRecord},
                                    p::TBLGParameters)
    n = length(sites)
    rows, cols, vals = Int[], Int[], Float64[]
    for i in 1:n, j in (i + 1):n
        t = hopping_value(sites[i], sites[j], p)
        iszero(t) && continue
        push!(rows, i); push!(cols, j); push!(vals, t)
        push!(rows, j); push!(cols, i); push!(vals, t)
    end
    return sparse(rows, cols, vals, n, n)
end

quantics_level(n::Integer) = n > 0 ? trailing_zeros(n) :
    throw(ArgumentError("site count must be positive"))

function validate_site_addressing(sites::AbstractVector{<:SiteRecord};
                                  require_power_of_two=true)
    n = length(sites)
    (!require_power_of_two || ispow2(n)) ||
        throw(ArgumentError("physical site count $n is not an exact power of two"))
    ids = [s.id for s in sites]
    sort(ids) == collect(1:n) || throw(ArgumentError("site ids are not a bijection"))
    keys = [(s.layer, s.sublattice, s.n1, s.n2) for s in sites]
    length(unique(keys)) == n || throw(ArgumentError("duplicate lattice site record"))
    return (site_count=n, qtt_bits=quantics_level(n), bijective=true)
end

function count_by_layer_sublattice(sites)
    counts = Dict{Tuple{Int,Int},Int}()
    for s in sites
        key = (s.layer, s.sublattice)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end

"""
Find an exact radius for a requested site count, if a shell of the discrete
flake reaches it. nothing means the count is skipped by a degenerate shell.
"""
function find_exact_radius(target::Integer, p::TBLGParameters=TBLGParameters();
                           radius_max::Real=200.0,
                           center::Tuple{<:Real,<:Real}=(0.0, 0.0))
    target > 0 || throw(ArgumentError("target must be positive"))
    c = (Float64(center[1]), Float64(center[2]))
    candidates = circular_flake(Float64(radius_max), p; center=c)
    radii = sort(unique(hypot(s.x - c[1], s.y - c[2]) for s in candidates))
    for radius in radii
        n_inside = Base.count(s -> hypot(s.x - c[1], s.y - c[2]) <=
                                      radius + 10eps(radius), candidates)
        n_inside == target && return radius
        n_inside > target && break
    end
    return nothing
end

end
