
export extract_hartree_mpo_1d
struct HartreeEvaluate1D
    sys::System
end

function (he::HartreeEvaluate1D)(i_float::Real)
    i = Int(i_float) +1

    v_hartree = 0.0
    L = length(he.sys.sites)

    neighbors = [i-1, i+1]

    for j in neighbors
        if i <= j <= 2^L
            n_j = MatrixChecker(he.sys.ρ, 
                                he.sys.sites,
                                 j,
                                 j,
                                 he.sys.bra_states,
                                 he.sys.ket_states
                                 )

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
    evaluator = HartreeEvaluate1D(sys)

    _, vh_mpo, _ = Quantics_TCI(
        x -> evaluator(x),
        Float64,
        sys.sites,
        sys.params.tci_tol
    )
    return vh_mpo
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
    ρ_val = MatrixChecker(fe.sys.ρ, fe.sys.sites, x+1, x, fe.sys.bra_states, fe.sys.ket_states)
    
    # The Fock coefficient is -U * bond_order
    return - fe.sys.params.U * real(ρ_val)
end


function extract_fock_mpo_1d(sys::System)
    sites = sys.sites
    params = sys.params
    

    fe = FockEvaluator1D(sys)
    _, F_MPO, _ = Quantics_TCI(x -> fe(x), Float64, sites, params.tci_tol)
    
    T_R, T_L = _build_translation_chain(sites) 
    

    VF_R = apply(F_MPO, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    VF_L = apply(T_L, ITensors.dag(F_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    
    return +(VF_R, VF_L; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
end


    
