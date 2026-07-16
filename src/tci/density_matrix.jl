
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


    
