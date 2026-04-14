
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
