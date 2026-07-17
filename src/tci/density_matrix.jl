
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
function _density_diagonal_qtt_tensors(rho::MPO, sites::Vector{<:Index})
    tensors = Vector{ITensor}(undef, length(sites))
    for site_number in eachindex(sites)
        site = sites[site_number]
        output_site = sim(site)
        projected = rho[site_number] * delta(prime(site), site, output_site)
        tensors[site_number] = replaceind(projected, output_site => site)
    end
    return tensors
end

function _diagonal_mpo_from_qtt_tensors(
    tensors::Vector{ITensor},
    sites::Vector{<:Index},
    params::AbstractModelParameters;
    symmetrize::Bool=true,
)
    diagonal = MPO(length(sites))
    for site_number in eachindex(sites)
        diagonal[site_number] = Quantics._asdiagonal(tensors[site_number], sites[site_number])
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
    extract_mean_fields(sys)

Return `(VH, VF)` using the geometry-specific nearest-neighbour Hartree/Fock
implementation for `sys`. Square `VF` is the sum of independently constructed
horizontal and vertical exchange fields.
"""
function _extract_mean_fields_with_components(sys::System)
    if sys.params isa Parameters1D
        return extract_hartree_mpo_1d(sys), extract_fock_mpo_1d(sys), nothing
    elseif sys.params isa ParametersSquare
        horizontal = extract_fock_mpo_square_horizontal(sys)
        vertical = extract_fock_mpo_square_vertical(sys)
        fock = +(horizontal, vertical;
            cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
        )
        return extract_hartree_mpo_square(sys), fock, (horizontal=horizontal, vertical=vertical)
    end
    throw(ArgumentError("no mean-field extractor for $(typeof(sys.params))"))
end

function extract_mean_fields(sys::System)
    hartree, fock, _ = _extract_mean_fields_with_components(sys)
    return hartree, fock
end


    
