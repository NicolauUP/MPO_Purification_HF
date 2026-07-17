
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
function extract_mean_fields(sys::System)
    if sys.params isa Parameters1D
        return extract_hartree_mpo_1d(sys), extract_fock_mpo_1d(sys)
    elseif sys.params isa ParametersSquare
        horizontal = extract_fock_mpo_square_horizontal(sys)
        vertical = extract_fock_mpo_square_vertical(sys)
        fock = +(horizontal, vertical;
            cutoff=sys.params.itensors_tol, maxdim=sys.params.itensors_maxdim,
        )
        return extract_hartree_mpo_square(sys), fock
    end
    throw(ArgumentError("no mean-field extractor for $(typeof(sys.params))"))
end


    
