
export extract_hartree_mpo_1d
struct HartreeEvaluate1D
    sys::System
    density_cache::Dict{Int,Float64}
end

function (he::HartreeEvaluate1D)(i_float::Real)
    i = Int(i_float) +1

    v_hartree = 0.0
    L = length(he.sys.sites)

    neighbors = [i-1, i+1]

    for j in neighbors
        if 1 <= j <= 2^L
            n_j = get!(he.density_cache, j) do
                real(MatrixChecker(
                    he.sys.ρ,
                    he.sys.sites,
                    j,
                    j,
                    he.sys.bra_states,
                    he.sys.ket_states,
                ))
            end

            v_hartree += n_j * he.sys.params.U 
        end
    end
       
    return v_hartree
end

"""
    extract_hartree_mpo_1d(sys::System)

Runs TCI to generate the 1D Hartree MPO.
"""

function extract_hartree_mpo_1d(sys::System)
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    evaluator = HartreeEvaluate1D(sys, Dict{Int,Float64}())
    return diagonal_mpo_from_function(
        x -> evaluator(x),
        Float64,
        sys.sites,
        sys.params.tci_tol
    )
end

"""
    _density_diagonal_mpo(rho, sites, params)

Project the bra/ket physical indices of an MPO onto equal local states and
return the resulting diagonal MPO. The symmetrization preserves the current
Hartree convention `real(rho[i,i])` when finite MPO truncation leaves a small
anti-Hermitian diagonal component.
"""
function _density_diagonal_qtt_tensors(
    rho::MPO,
    sites::Vector{<:Index};
    to_backend=identity,
)
    tensors = Vector{ITensor}(undef, length(sites))
    for site_number in eachindex(sites)
        site = sites[site_number]
        output_site = sim(site)
        projected = rho[site_number] * to_backend(delta(prime(site), site, output_site))
        tensors[site_number] = replaceind(projected, output_site => site)
    end
    return tensors
end

function _diagonal_mpo_from_qtt_tensors(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    params::AbstractModelParameters;
    symmetrize::Bool=true,
    to_backend=identity,
)
    diagonal = MPO(length(sites))
    for site_number in eachindex(sites)
        if to_backend === identity
            diagonal[site_number] = Quantics._asdiagonal(
                tensors[site_number], sites[site_number],
            )
        else
            # Quantics._asdiagonal currently materializes `Array(tensor)`,
            # which forces device tensors through the host. The three-index
            # delta duplicates the coefficient's physical bit into the MPO
            # bra/ket indices while preserving the tensor storage backend.
            site = sites[site_number]
            output_site = sim(site)
            diagonalizer = to_backend(delta(site, prime(site), output_site))
            diagonal[site_number] = replaceind(
                tensors[site_number] * diagonalizer,
                output_site => site,
            )
        end
    end
    if symmetrize
        return +(0.5 * diagonal, 0.5 * ITensors.dag(diagonal);
            cutoff=params.itensors_tol,
            maxdim=params.itensors_maxdim,
        )
    end
    return diagonal
end

function _density_diagonal_mpo(
    rho::MPO,
    sites::Vector{<:Index},
    params::AbstractModelParameters,
)
    tensors = _density_diagonal_qtt_tensors(rho, sites)
    return _diagonal_mpo_from_qtt_tensors(tensors, sites, params)
end

function _shift_diagonal_mpo(
    diagonal::MPO,
    forward::MPO,
    backward::MPO,
    params::AbstractModelParameters,
)
    shifted = apply(forward, diagonal;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return apply(shifted, backward;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
end

function _binary_shift_transition(output_bit::Int, carry_in::Int, direction::Symbol)
    if direction == :right
        value = output_bit + carry_in
        return value % 2, value ÷ 2
    elseif direction == :left
        value = output_bit - carry_in
        return mod(value, 2), value < 0 ? 1 : 0
    end
    throw(ArgumentError("direction must be :right or :left, got $direction"))
end

"""
    _shift_qtt_tensors_binary_carry(tensors, sites, direction)

Apply the open-boundary binary shift `f(i) -> f(i+1)` (`:right`) or
`f(i) -> f(i-1)` (`:left`) directly to QTT/MPS tensors. The carry enters at
the least-significant bit and is fused into the existing QTT bond at every
internal cut. No general translation-MPO product is formed.
"""
function _shift_qtt_tensors_binary_carry(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    direction::Symbol,
)
    length(tensors) == length(sites) || throw(ArgumentError(
        "tensor and site counts must agree",
    ))
    L = length(sites)
    carries = [Index(2, "HartreeCarry,link=$link") for link in 1:(L - 1)]
    shifted = Vector{ITensor}(undef, L)

    for site_number in 1:L
        local_tensor = nothing
        for output_bit in 0:1, carry_in in 0:1
            input_bit, carry_out = _binary_shift_transition(
                output_bit, carry_in, direction,
            )
            site_number == L && carry_in != 1 && continue
            site_number == 1 && carry_out != 0 && continue

            term = tensors[site_number] * onehot(sites[site_number] => input_bit + 1)
            term *= onehot(sites[site_number] => output_bit + 1)
            site_number > 1 && (term *= onehot(carries[site_number - 1] => carry_out + 1))
            site_number < L && (term *= onehot(carries[site_number] => carry_in + 1))
            local_tensor = isnothing(local_tensor) ? term : local_tensor + term
        end
        shifted[site_number] = local_tensor
    end

    for site_number in 1:(L - 1)
        density_link = commonind(tensors[site_number], tensors[site_number + 1])
        carry_link = carries[site_number]
        fuse_links = combiner(density_link, carry_link)
        shifted[site_number] *= fuse_links
        shifted[site_number + 1] *= ITensors.dag(fuse_links)
    end
    return shifted
end

"""
    _density_superdiagonal_qtt_tensors_binary_carry(rho, sites)

Contract the density MPO locally onto its open-chain superdiagonal,
`rho[i, i + 1]`. The QTT output coordinate is the row index `i`; a binary
increment carry selects the column index `i + 1`. The carry enters at the
least-significant bit, and an overflow out of the most-significant bit is
discarded, leaving the final output coefficient exactly zero.
"""
function _density_superdiagonal_qtt_tensors_binary_carry(
    rho::MPO,
    sites::Vector{<:Index},
)
    length(rho) == length(sites) || throw(ArgumentError(
        "density MPO and site counts must agree",
    ))
    L = length(sites)
    carries = [Index(2, "FockCarry,link=$link") for link in 1:(L - 1)]
    band = Vector{ITensor}(undef, L)

    for site_number in 1:L
        site = sites[site_number]
        output_site = sim(site)
        local_tensor = nothing
        for row_bit in 0:1, carry_in in 0:1
            column_bit, carry_out = _binary_shift_transition(
                row_bit, carry_in, :right,
            )
            site_number == L && carry_in != 1 && continue
            site_number == 1 && carry_out != 0 && continue

            term = rho[site_number] * onehot(prime(site) => row_bit + 1)
            term *= onehot(site => column_bit + 1)
            term *= onehot(output_site => row_bit + 1)
            site_number > 1 && (term *= onehot(carries[site_number - 1] => carry_out + 1))
            site_number < L && (term *= onehot(carries[site_number] => carry_in + 1))
            local_tensor = isnothing(local_tensor) ? term : local_tensor + term
        end
        band[site_number] = replaceind(local_tensor, output_site => site)
    end

    for site_number in 1:(L - 1)
        density_link = commonind(rho[site_number], rho[site_number + 1])
        carry_link = carries[site_number]
        fuse_links = combiner(density_link, carry_link)
        band[site_number] *= fuse_links
        band[site_number + 1] *= ITensors.dag(fuse_links)
    end
    return band
end

"""
    extract_hartree_mpo_binary_carry_1d(sys)

Experimental fused tensorial 1D Hartree extractor. It obtains the QTT density
diagonal by local bra/ket contraction and applies binary increment/decrement
carry transitions directly to those tensors. It avoids scalar TCI queries and
the general MPO products used by `extract_hartree_mpo_tensorial_1d`.
"""
function extract_hartree_mpo_binary_carry_1d(sys::System)
    sys.params isa Parameters1D || throw(ArgumentError(
        "extract_hartree_mpo_binary_carry_1d requires Parameters1D",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    density_tensors = _density_diagonal_qtt_tensors(sys.ρ, sys.sites)
    right_tensors = _shift_qtt_tensors_binary_carry(density_tensors, sys.sites, :right)
    left_tensors = _shift_qtt_tensors_binary_carry(density_tensors, sys.sites, :left)
    right_density = _diagonal_mpo_from_qtt_tensors(
        right_tensors, sys.sites, sys.params; symmetrize=false,
    )
    left_density = _diagonal_mpo_from_qtt_tensors(
        left_tensors, sys.sites, sys.params; symmetrize=false,
    )
    hartree = +(sys.params.U * right_density, sys.params.U * left_density;
        cutoff=sys.params.itensors_tol,
        maxdim=sys.params.itensors_maxdim,
    )
    return +(0.5 * hartree, 0.5 * ITensors.dag(hartree);
        cutoff=sys.params.itensors_tol,
        maxdim=sys.params.itensors_maxdim,
    )
end

"""
    _shift_qtt_tensors_binary_carry_square(tensors, sites, direction)

Shift a diagonal QTT field by one open-square neighbour without forming a
translation MPO. Square coordinates use the established interleaved bit order
`(y₀, x₀, y₁, x₁, …)`. QTT tensor order runs from the most-significant bit to
the least-significant bit, so horizontal shifts carry through odd Julia site
positions (`1, 3, …`) and vertical shifts through even positions (`2, 4, …`).
The other coordinate's bits transmit the carry state unchanged. Overflow and
underflow are discarded at the physical square boundary.
"""
function _shift_qtt_tensors_binary_carry_square(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    direction::Symbol,
    ;
    to_backend=identity,
)
    length(tensors) == length(sites) || throw(ArgumentError(
        "tensor and site counts must agree",
    ))
    L = length(sites)
    iseven(L) || throw(ArgumentError("square QTT shift requires an even number of sites"))
    direction in (:right, :left, :up, :down) || throw(ArgumentError(
        "square direction must be :right, :left, :up, or :down; got $direction",
    ))

    horizontal = direction in (:right, :left)
    binary_direction = direction in (:right, :up) ? :right : :left
    active_site(site_number) = horizontal ? isodd(site_number) : iseven(site_number)
    carries = [Index(2, "SquareHartreeCarry,link=$link") for link in 1:(L - 1)]
    shifted = Vector{ITensor}(undef, L)

    for site_number in 1:L
        site = sites[site_number]
        local_tensor = nothing
        is_active = active_site(site_number)
        for output_bit in 0:1, carry_in in 0:1
            input_bit, carry_out = if is_active
                _binary_shift_transition(output_bit, carry_in, binary_direction)
            else
                output_bit, carry_in
            end

            # The carry begins at the most-significant bit of the selected
            # coordinate. When the final QTT site belongs to the other axis,
            # it injects that initial carry into the next (active) site.
            if site_number == L
                is_active ? carry_in == 1 || continue : carry_out == 1 || continue
            end
            # The carry must vanish below the least-significant selected bit.
            # For a horizontal shift that bit is site 2, so site 1 consumes
            # the final zero carry without changing its (y) bit.
            if site_number == 1
                is_active ? carry_out == 0 || continue : carry_in == 0 || continue
            end

            term = tensors[site_number] * to_backend(onehot(site => input_bit + 1))
            term *= to_backend(onehot(site => output_bit + 1))
            site_number > 1 && (term *= to_backend(
                onehot(carries[site_number - 1] => carry_out + 1),
            ))
            site_number < L && (term *= to_backend(
                onehot(carries[site_number] => carry_in + 1),
            ))
            local_tensor = isnothing(local_tensor) ? term : local_tensor + term
        end
        shifted[site_number] = local_tensor
    end

    for site_number in 1:(L - 1)
        density_link = commonind(tensors[site_number], tensors[site_number + 1])
        fuse_links = to_backend(combiner(density_link, carries[site_number]))
        shifted[site_number] *= fuse_links
        shifted[site_number + 1] *= ITensors.dag(fuse_links)
    end
    return shifted
end

"""
    _density_neighbour_band_qtt_tensors_binary_carry_square(rho, sites, direction)

Locally contract an interleaved square density MPO onto a forward nearest-
neighbour band. `direction=:right` returns the QTT coefficient field
`rho[(x,y), (x+1,y)]`; `direction=:up` returns
`rho[(x,y), (x,y+1)]`. Open-boundary bonds have exactly zero coefficients.

The carry acts only on the selected coordinate's QTT bits: odd Julia tensor
positions for `:right` and even positions for `:up`. It is the square analogue
of `_density_superdiagonal_qtt_tensors_binary_carry`, preserving the density
matrix's row/column orientation until the caller applies the established
real-exchange convention.
"""
function _density_neighbour_band_qtt_tensors_binary_carry_square(
    rho::MPO,
    sites::Vector{<:Index},
    direction::Symbol,
    ;
    to_backend=identity,
)
    length(rho) == length(sites) || throw(ArgumentError(
        "density MPO and site count must agree",
    ))
    L = length(sites)
    iseven(L) || throw(ArgumentError("square QTT band extraction requires an even number of sites"))
    direction in (:right, :up) || throw(ArgumentError(
        "square Fock band direction must be :right or :up; got $direction",
    ))

    horizontal = direction == :right
    active_site(site_number) = horizontal ? isodd(site_number) : iseven(site_number)
    carries = [Index(2, "SquareFockCarry,link=$link") for link in 1:(L - 1)]
    band = Vector{ITensor}(undef, L)

    for site_number in 1:L
        site = sites[site_number]
        output_site = sim(site)
        local_tensor = nothing
        is_active = active_site(site_number)
        for row_bit in 0:1, carry_in in 0:1
            column_bit, carry_out = if is_active
                _binary_shift_transition(row_bit, carry_in, :right)
            else
                row_bit, carry_in
            end

            # Inject the unit increment at the most-significant selected bit,
            # then require it to vanish below the selected coordinate's least
            # significant bit. The other axis only transports the carry.
            if site_number == L
                is_active ? carry_in == 1 || continue : carry_out == 1 || continue
            end
            if site_number == 1
                is_active ? carry_out == 0 || continue : carry_in == 0 || continue
            end

            term = rho[site_number] * to_backend(onehot(prime(site) => row_bit + 1))
            term *= to_backend(onehot(site => column_bit + 1))
            term *= to_backend(onehot(output_site => row_bit + 1))
            site_number > 1 && (term *= to_backend(
                onehot(carries[site_number - 1] => carry_out + 1),
            ))
            site_number < L && (term *= to_backend(
                onehot(carries[site_number] => carry_in + 1),
            ))
            local_tensor = isnothing(local_tensor) ? term : local_tensor + term
        end
        band[site_number] = replaceind(local_tensor, output_site => site)
    end

    for site_number in 1:(L - 1)
        density_link = commonind(rho[site_number], rho[site_number + 1])
        fuse_links = to_backend(combiner(density_link, carries[site_number]))
        band[site_number] *= fuse_links
        band[site_number + 1] *= ITensors.dag(fuse_links)
    end
    return band
end

"""
    square_neighbour_adjacency_mpo(sites)

Return the exact QTT/MPO representation of the open-square nearest-neighbour
adjacency operator

```math
A_{\\mathrm{nn}} = S_x^+ + S_x^- + S_y^+ + S_y^- .
```

The virtual state is `(direction, carry)`: four direction branches and the
binary carry/borrow required by each branch. This produces the four-neighbour
sum in one MPO--MPS application, rather than forming four shifted fields and
compressing their MPO sum. Tensor order is most- to least-significant QTT bit,
with the established interleaved square convention `(y₀, x₀, y₁, x₁, ...)`.
"""
function square_neighbour_adjacency_mpo(sites::Vector{<:Index})
    L = length(sites)
    iseven(L) || throw(ArgumentError(
        "square QTT adjacency requires an even number of sites",
    ))
    directions = (:right, :left, :up, :down)
    state_links = [Index(8, "SquareAdjacencyState,link=$link") for link in 1:(L - 1)]
    adjacency = MPO(L)

    for site_number in 1:L
        site = sites[site_number]
        output_site = prime(site)
        local_tensor = nothing
        for (direction_number, direction) in enumerate(directions)
            horizontal = direction in (:right, :left)
            binary_direction = direction in (:right, :up) ? :right : :left
            active = horizontal ? isodd(site_number) : iseven(site_number)
            for output_bit in 0:1, carry_in in 0:1
                input_bit, carry_out = active ? _binary_shift_transition(
                    output_bit, carry_in, binary_direction,
                ) : (output_bit, carry_in)

                # At the least-significant physical coordinate bit, seed the
                # increment/decrement with one. At the most-significant bit,
                # retain only paths whose carry has vanished, which exactly
                # enforces the open boundary.
                site_number == L && carry_in != 1 && continue
                site_number == 1 && carry_out != 0 && continue

                term = onehot(output_site => output_bit + 1)
                term *= onehot(site => input_bit + 1)
                site_number > 1 && (term *= onehot(
                    state_links[site_number - 1] => 2 * (direction_number - 1) + carry_out + 1,
                ))
                site_number < L && (term *= onehot(
                    state_links[site_number] => 2 * (direction_number - 1) + carry_in + 1,
                ))
                local_tensor = isnothing(local_tensor) ? term : local_tensor + term
            end
        end
        adjacency[site_number] = local_tensor
    end
    return adjacency
end

"""
    _apply_square_neighbour_adjacency(tensors, sites; cutoff, maxdim)

Apply [`square_neighbour_adjacency_mpo`](@ref) to the density QTT. The
compression policy is explicit because the raw product can have a much larger
bond dimension than the physical Hartree field.
"""
function _apply_square_neighbour_adjacency(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    ;
    cutoff::Real,
    maxdim::Integer,
    to_backend=identity,
)
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError(
        "square Hartree cutoff must be finite and nonnegative",
    ))
    maxdim > 0 || throw(ArgumentError("square Hartree maxdim must be positive"))
    density = MPS(tensors)
    adjacency = to_backend(square_neighbour_adjacency_mpo(sites))
    field = apply(adjacency, density;
        cutoff=cutoff,
        maxdim=maxdim,
    )
    return copy(field.data)
end

"""
    extract_hartree_mpo_binary_carry_square(sys)

Construct the open-square nearest-neighbour Hartree field by projecting the
density diagonal once and applying exact binary-carry shifts along interleaved
`x` and `y` QTT bits. It avoids generic translation-MPO products, so a missing
physical neighbour is represented by a discarded carry overflow/underflow
rather than a truncation-sensitive boundary mask.
"""
function extract_hartree_mpo_binary_carry_square(sys::System)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_hartree_mpo_binary_carry_square requires ParametersSquare",
    ))
    params = sys.params
    iszero(params.U) && return zero_mpo(sys.sites)
    density_tensors = _density_diagonal_qtt_tensors(sys.ρ, sys.sites)
    shifted = (
        _shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :right),
        _shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :left),
        _shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :up),
        _shift_qtt_tensors_binary_carry_square(density_tensors, sys.sites, :down),
    )
    fields = map(tensors -> _diagonal_mpo_from_qtt_tensors(
        tensors, sys.sites, params; symmetrize=false,
    ), shifted)
    hartree = +(map(field -> params.U * field, fields)...;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return +(0.5 * hartree, 0.5 * ITensors.dag(hartree);
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
end

"""
    extract_hartree_mpo_binary_carry_square_adjacency(
        sys; cutoff=1e-12, maxdim=64,
    )

Fused square Hartree extractor. It constructs the density diagonal once and
evaluates `U * A_nn * n` in a single QTT MPO--MPS application.
Unlike [`extract_hartree_mpo_binary_carry_square`](@ref), it does not form or
compress a generic sum of four global diagonal MPOs.

The default compression policy is based on the million-site square benchmark:
`cutoff=1e-12, maxdim=64` reduced the raw field from bond dimension 512 to 17
while keeping the largest sampled local-field error below `5e-5`. Pass
`cutoff=0, maxdim=typemax(Int)` to recover the untruncated diagnostic path.
This is the square Hartree path used by [`extract_mean_fields`](@ref) and
`run_scf!`.
"""
function extract_hartree_mpo_binary_carry_square_adjacency(
    sys::System;
    cutoff::Real=1e-12,
    maxdim::Integer=64,
    to_backend=identity,
)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_hartree_mpo_binary_carry_square_adjacency requires ParametersSquare",
    ))
    params = sys.params
    iszero(params.U) && return zero_mpo(sys.sites)
    density_tensors = _density_diagonal_qtt_tensors(
        sys.ρ, sys.sites; to_backend=to_backend,
    )
    field_tensors = _apply_square_neighbour_adjacency(
        density_tensors, sys.sites;
        cutoff=cutoff, maxdim=maxdim, to_backend=to_backend,
    )
    field = _diagonal_mpo_from_qtt_tensors(
        field_tensors, sys.sites, params;
        symmetrize=false, to_backend=to_backend,
    )
    return params.U * field
end

"""
    extract_fock_mpo_binary_carry_square_horizontal(sys)

Square horizontal Fock extractor. It obtains the right-directed
band `rho[(x,y),(x+1,y)]` by a local interleaved binary carry, applies the
current `-U * real(rho[i,j])` convention, and assembles the established
Hermitian horizontal field. This is the default square SCF path; the cached
TCI extractor remains available through `square_fock_method=:tci`.
"""
function extract_fock_mpo_binary_carry_square_horizontal(
    sys::System;
    to_backend=identity,
    translations=sys.translations,
)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_fock_mpo_binary_carry_square_horizontal requires ParametersSquare",
    ))
    params = sys.params
    iszero(params.U) && return zero_mpo(sys.sites)
    band_tensors = _density_neighbour_band_qtt_tensors_binary_carry_square(
        sys.ρ, sys.sites, :right; to_backend=to_backend,
    )
    coefficients = -params.U * _diagonal_mpo_from_qtt_tensors(
        band_tensors, sys.sites, params;
        symmetrize=true, to_backend=to_backend,
    )
    T_R, T_L, _, _ = translations
    right_term = apply(coefficients, T_R;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    left_term = apply(T_L, ITensors.dag(coefficients);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    return +(right_term, left_term;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

"""
    extract_fock_mpo_binary_carry_square_vertical(sys)

Square vertical counterpart of
[`extract_fock_mpo_binary_carry_square_horizontal`](@ref). It samples only
the up-directed band `rho[(x,y),(x,y+1)]` and applies the same real-exchange
and open-boundary conventions.
"""
function extract_fock_mpo_binary_carry_square_vertical(
    sys::System;
    to_backend=identity,
    translations=sys.translations,
)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_fock_mpo_binary_carry_square_vertical requires ParametersSquare",
    ))
    params = sys.params
    iszero(params.U) && return zero_mpo(sys.sites)
    band_tensors = _density_neighbour_band_qtt_tensors_binary_carry_square(
        sys.ρ, sys.sites, :up; to_backend=to_backend,
    )
    coefficients = -params.U * _diagonal_mpo_from_qtt_tensors(
        band_tensors, sys.sites, params;
        symmetrize=true, to_backend=to_backend,
    )
    _, _, T_U, T_D = translations
    up_term = apply(coefficients, T_U;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    down_term = apply(T_D, ITensors.dag(coefficients);
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
    return +(up_term, down_term;
        cutoff=params.itensors_tol, maxdim=params.itensors_maxdim,
    )
end

"""
    extract_fock_mpo_binary_carry_1d(sys)

Experimental 1D Fock extractor. It contracts `rho[i, i + 1]` directly with a
local binary carry, converts the result to the established real exchange
coefficient `-U * real(rho[i, i + 1])`, then uses the existing translation-MPO
assembly for the Hermitian nearest-neighbour field. Cached TCI remains the
default path.
"""
function extract_fock_mpo_binary_carry_1d(sys::System)
    sys.params isa Parameters1D || throw(ArgumentError(
        "extract_fock_mpo_binary_carry_1d requires Parameters1D",
    ))
    params = sys.params
    iszero(params.U) && return zero_mpo(sys.sites)

    band_tensors = _density_superdiagonal_qtt_tensors_binary_carry(sys.ρ, sys.sites)
    real_bond_order = _diagonal_mpo_from_qtt_tensors(
        band_tensors, sys.sites, params; symmetrize=true,
    )
    coefficients = -params.U * real_bond_order
    T_R, T_L = sys.translations
    right_term = apply(coefficients, T_R;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    left_term = apply(T_L, ITensors.dag(coefficients);
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
    return +(right_term, left_term;
        cutoff=params.itensors_tol,
        maxdim=params.itensors_maxdim,
    )
end

"""
    extract_hartree_mpo_tensorial_1d(sys)

Experimental tensorial 1D Hartree extractor. It projects the density MPO
diagonal directly, then shifts it with the open-chain translation MPOs:
`U * (T_R D_n T_L + T_L D_n T_R)`. This avoids scalar `MatrixChecker` calls
and TCI reconstruction, but remains opt-in until benchmarks establish its
runtime and bond-dimension behavior.
"""
function extract_hartree_mpo_tensorial_1d(sys::System)
    sys.params isa Parameters1D || throw(ArgumentError(
        "extract_hartree_mpo_tensorial_1d requires Parameters1D",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    density_diagonal = _density_diagonal_mpo(sys.ρ, sys.sites, sys.params)
    T_R, T_L = sys.translations
    right_density = _shift_diagonal_mpo(density_diagonal, T_R, T_L, sys.params)
    left_density = _shift_diagonal_mpo(density_diagonal, T_L, T_R, sys.params)
    return +(sys.params.U * right_density, sys.params.U * left_density;
        cutoff=sys.params.itensors_tol,
        maxdim=sys.params.itensors_maxdim,
    )
end

struct HartreeEvaluateSquare
    sys::System
    density_cache::Dict{Int,Float64}
end

function (evaluator::HartreeEvaluateSquare)(coordinate::Real)
    site = Int(coordinate) + 1
    value = 0.0
    for neighbour in values(square_neighbours(site, evaluator.sys.params.L))
        isnothing(neighbour) && continue
        density = get!(evaluator.density_cache, neighbour) do
            real(MatrixChecker(
                evaluator.sys.ρ,
                evaluator.sys.sites,
                neighbour,
                neighbour,
                evaluator.sys.bra_states,
                evaluator.sys.ket_states,
            ))
        end
        value += evaluator.sys.params.U * density
    end
    return value
end

"""
    extract_hartree_mpo_square(sys)

Construct the open-boundary square-lattice Hartree field. At each site it is
`U` times the sum of the density on its valid right, left, up, and down
neighbours. No periodic wraparound terms are included.
"""
function extract_hartree_mpo_square(sys::System)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_hartree_mpo_square requires ParametersSquare",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    evaluator = HartreeEvaluateSquare(sys, Dict{Int,Float64}())
    return diagonal_mpo_from_function(
        x -> evaluator(x), Float64, sys.sites, sys.params.tci_tol,
    )
end

"""
    extract_hartree_mpo_tensorial_square(sys)

Experimental tensorial square Hartree extractor. It directly projects the
density diagonal and applies the four open-boundary translations, preserving
the same four-neighbour convention as [`extract_hartree_mpo_square`](@ref).
"""
function extract_hartree_mpo_tensorial_square(sys::System)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_hartree_mpo_tensorial_square requires ParametersSquare",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    density_diagonal = _density_diagonal_mpo(sys.ρ, sys.sites, sys.params)
    T_R, T_L, T_U, T_D = sys.translations
    right_density = _shift_diagonal_mpo(density_diagonal, T_R, T_L, sys.params)
    left_density = _shift_diagonal_mpo(density_diagonal, T_L, T_R, sys.params)
    up_density = _shift_diagonal_mpo(density_diagonal, T_U, T_D, sys.params)
    down_density = _shift_diagonal_mpo(density_diagonal, T_D, T_U, sys.params)
    return +(sys.params.U * right_density, sys.params.U * left_density,
        sys.params.U * up_density, sys.params.U * down_density;
        cutoff=sys.params.itensors_tol,
        maxdim=sys.params.itensors_maxdim,
    )
end



struct FockEvaluator1D
    sys::System
end

function (fe::FockEvaluator1D)(i_float::Real)
    # i_float is the coordinate x
    x = Int(i_float) + 1
    L = fe.sys.params.L
    
    # Boundary: No bond exists to the right of the last site
    if x >= 2^L
        return 0.0
    end
    
    # Measure the bond order <c†_{x+1} c_x>
    # We use your MatrixChecker logic
    ρ_val = MatrixChecker(fe.sys.ρ, fe.sys.sites, x, x+1, fe.sys.bra_states, fe.sys.ket_states)
    
    # The Fock coefficient is -U * bond_order
    return - fe.sys.params.U * real(ρ_val)
end


function extract_fock_mpo_1d(sys::System)
    sites = sys.sites
    params = sys.params
    

    iszero(params.U) && return zero_mpo(sites)
    fe = FockEvaluator1D(sys)
    F_MPO = diagonal_mpo_from_function(x -> fe(x), Float64, sites, params.tci_tol)
    
    T_R, T_L = sys.translations
    

    VF_R = apply(F_MPO, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    VF_L = apply(T_L, ITensors.dag(F_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    
    return +(VF_R, VF_L; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
end

struct HorizontalFockEvaluatorSquare
    sys::System
end

function (evaluator::HorizontalFockEvaluatorSquare)(coordinate::Real)
    site = Int(coordinate) + 1
    right = square_neighbours(site, evaluator.sys.params.L).right
    isnothing(right) && return 0.0
    bond_order = MatrixChecker(
        evaluator.sys.ρ,
        evaluator.sys.sites,
        site,
        right,
        evaluator.sys.bra_states,
        evaluator.sys.ket_states,
    )
    return -evaluator.sys.params.U * real(bond_order)
end

"""
    extract_fock_mpo_square_horizontal(sys)

Construct the real-exchange Fock field on right-directed open-square bonds.
For every horizontal bond `i -> j`, its coefficient is
`-U * real(rho[i,j])`, and the returned MPO includes the Hermitian conjugate
left-directed term. Vertical exchange is intentionally not included here.
"""
function extract_fock_mpo_square_horizontal(sys::System)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_fock_mpo_square_horizontal requires ParametersSquare",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    evaluator = HorizontalFockEvaluatorSquare(sys)
    coefficients = diagonal_mpo_from_function(
        x -> evaluator(x), Float64, sys.sites, sys.params.tci_tol,
    )
    T_R, T_L, _, _ = sys.translations
    right_term = apply(coefficients, T_R;
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
    left_term = apply(T_L, ITensors.dag(coefficients);
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
    return +(right_term, left_term;
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
end

struct VerticalFockEvaluatorSquare
    sys::System
end

function (evaluator::VerticalFockEvaluatorSquare)(coordinate::Real)
    site = Int(coordinate) + 1
    up = square_neighbours(site, evaluator.sys.params.L).up
    isnothing(up) && return 0.0
    bond_order = MatrixChecker(
        evaluator.sys.ρ,
        evaluator.sys.sites,
        site,
        up,
        evaluator.sys.bra_states,
        evaluator.sys.ket_states,
    )
    return -evaluator.sys.params.U * real(bond_order)
end

"""
    extract_fock_mpo_square_vertical(sys)

Construct the real-exchange Fock field on up-directed open-square bonds. For
every vertical bond `i -> j`, its coefficient is `-U * real(rho[i,j])`, and
the returned MPO includes the Hermitian conjugate down-directed term.
Horizontal exchange is intentionally not included here.
"""
function extract_fock_mpo_square_vertical(sys::System)
    sys.params isa ParametersSquare || throw(ArgumentError(
        "extract_fock_mpo_square_vertical requires ParametersSquare",
    ))
    iszero(sys.params.U) && return zero_mpo(sys.sites)
    evaluator = VerticalFockEvaluatorSquare(sys)
    coefficients = diagonal_mpo_from_function(
        x -> evaluator(x), Float64, sys.sites, sys.params.tci_tol,
    )
    _, _, T_U, T_D = sys.translations
    up_term = apply(coefficients, T_U;
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
    down_term = apply(T_D, ITensors.dag(coefficients);
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
    return +(up_term, down_term;
        cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
    )
end

"""
    extract_mean_fields(sys; square_fock_method=:binary_carry)

Return `(VH, VF)` using the geometry-specific nearest-neighbour Hartree/Fock
implementation for `sys`. Square `VF` is the sum of independently constructed
horizontal and vertical exchange fields. `square_fock_method` may be `:tci`
or `:binary_carry`; it has no effect for one-dimensional systems.
"""
function _extract_mean_fields_with_components(
    sys::System;
    square_fock_method::Symbol=:binary_carry,
    phase_callback::Union{Nothing,Function}=nothing,
    phase_synchronize::Function=() -> nothing,
)
    square_fock_method in (:tci, :binary_carry) || throw(ArgumentError(
        "square_fock_method must be :tci or :binary_carry, got $square_fock_method",
    ))
    if sys.params isa Parameters1D
        hartree = _profiled_operation(
            :hartree; callback=phase_callback, synchronize=phase_synchronize,
        ) do
            extract_hartree_mpo_binary_carry_1d(sys)
        end
        fock = _profiled_operation(
            :fock; callback=phase_callback, synchronize=phase_synchronize,
        ) do
            extract_fock_mpo_binary_carry_1d(sys)
        end
        return hartree, fock, nothing
    elseif sys.params isa ParametersSquare
        if square_fock_method == :tci
            horizontal = _profiled_operation(
                :fock_horizontal; callback=phase_callback, synchronize=phase_synchronize,
            ) do
                extract_fock_mpo_square_horizontal(sys)
            end
            vertical = _profiled_operation(
                :fock_vertical; callback=phase_callback, synchronize=phase_synchronize,
            ) do
                extract_fock_mpo_square_vertical(sys)
            end
        else
            horizontal = _profiled_operation(
                :fock_horizontal; callback=phase_callback, synchronize=phase_synchronize,
            ) do
                extract_fock_mpo_binary_carry_square_horizontal(sys)
            end
            vertical = _profiled_operation(
                :fock_vertical; callback=phase_callback, synchronize=phase_synchronize,
            ) do
                extract_fock_mpo_binary_carry_square_vertical(sys)
            end
        end
        fock = _profiled_operation(
            :fock_sum; callback=phase_callback, synchronize=phase_synchronize,
        ) do
            +(horizontal, vertical;
                cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
            )
        end
        hartree = _profiled_operation(
            :hartree; callback=phase_callback, synchronize=phase_synchronize,
        ) do
            extract_hartree_mpo_binary_carry_square_adjacency(sys)
        end
        return hartree, fock, (horizontal=horizontal, vertical=vertical)
    end
    throw(ArgumentError("no mean-field extractor for $(typeof(sys.params))"))
end

function extract_mean_fields(sys::System; square_fock_method::Symbol=:binary_carry)
    hartree, fock, _ = _extract_mean_fields_with_components(
        sys; square_fock_method=square_fock_method,
    )
    return hartree, fock
end


    
