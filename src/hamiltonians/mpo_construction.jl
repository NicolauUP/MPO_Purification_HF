# src/hamiltonian/mpo_construction.jl



function build_translation_chain(sites)
    L = length(sites)
    T_R_opsum = OpSum()
    T_L_opsum = OpSum()

    for l in 1:L
        opsum_R_temp = OpSum()
        opsum_L_temp = OpSum()

        opsum_R_temp += 1, "σ+", l
        opsum_L_temp += 1, "σ-", l

        for m in l+1:L
            opsum_R_temp *= 1, "σ-", m
            opsum_L_temp *= 1, "σ+", m

        end

        T_R_opsum += opsum_R_temp
        T_L_opsum += opsum_L_temp

    end

    T_R_MPO = MPO(T_R_opsum, sites)
    T_L_MPO = MPO(T_L_opsum, sites)

    return T_R_MPO, T_L_MPO
end

function build_translation_square(sites)
    L = length(sites)
    Lx = div(L, 2)
    Ly = div(L, 2)

    T_R_opsum = OpSum() #Hopping to the right
    T_L_opsum = OpSum() #Hopping to the left
    T_U_opsum = OpSum() #Hopping up
    T_D_opsum = OpSum() #Hopping down
    

    for l in 1:2:L
        opsum_R_temp = OpSum()
        opsum_L_temp = OpSum()
        opsum_U_temp = OpSum()
        opsum_D_temp = OpSum()

        opsum_R_temp += 1, "σ+", l
        opsum_L_temp += 1, "σ-", l
        opsum_U_temp += 1, "σ+", l+1
        opsum_D_temp += 1, "σ-", l+1
        

        for m in l+2:2:L
            opsum_R_temp *= 1, "σ-", m
            opsum_L_temp *= 1, "σ+", m
            opsum_U_temp *= 1, "σ-", m+1
            opsum_D_temp *= 1, "σ+", m+1
        end

        T_R_opsum += opsum_R_temp
        T_L_opsum += opsum_L_temp
        T_U_opsum += opsum_U_temp
        T_D_opsum += opsum_D_temp


    end
    T_R_MPO = MPO(T_R_opsum, sites)
    T_L_MPO = MPO(T_L_opsum, sites)
    T_U_MPO = MPO(T_U_opsum, sites)
    T_D_MPO = MPO(T_D_opsum, sites)

    



end



function build_W(sites, params)
    if !isnothing(params.W)
        _, W_MPO, _ = Quantics_TCI(params.W, Float64, sites, params.tci_tol)
        return W_MPO
    else
        return nothing
    end
end




function build_H0(sites, params::Parameters1D)
    T_R, T_L = build_translation_chain(sites)
 

    H0 = nothing
    if params.t isa Number
        println("Using constant hopping t = $(params.t)")

        H0 = params.t * (T_R + T_L)

    elseif params.t isa Function

        _, T_MPO, _ = Quantics_TCI(params.t, Float64, sites, params.tci_tol)

        H_T_R = apply(T_MPO, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        H_T_L = apply(T_L, ITensors.dag(T_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        H0 = +(H_T_R, H_T_L; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        
    end


    if !isnothing(params.W)
        W_MPO = build_W(sites, params)
        if isnothing(H0)
            H0 = W_MPO
        else
            H0 = +(H0, W_MPO; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        end
    end
    return H0
end



function build_H0(sites, params::ParametersSquare)
    println("Building 2D Hamiltonian MPO...")
    T_R, T_L, T_U, T_D = build_translation_square(sites)
 

    H0 = nothing
    tx, ty = params.t
    if tx isa Number && ty isa Number #They will always be the same type!
        println("Using constant hopping tx = $(tx) and ty = $(ty)")

        H0 = +(tx * (T_R + T_L), ty * (T_U + T_D); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

        # H0 = +(H0, ty * (T_U + T_D); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    

    elseif tx isa Function && ty isa Function

        _, TxMPO, _ = Quantics_TCI(tx, Float64, sites, params.tci_tol)
        _, TyMPO, _ = Quantics_TCI(ty, Float64, sites, params.tci_tol)

        H_T_R = apply(Tx_MPO, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        H_T_L = apply(T_L, ITensors.dag(Tx_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        H_T_U = apply(Ty_MPO, T_U; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        H_T_D = apply(T_D, ITensors.dag(Ty_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)

        H0 = +(H_T_R, H_T_L, H_T_U, H_T_D; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    end
        


    if !isnothing(params.W)
        W_MPO = build_W(sites, params)
        if isnothing(H0)
            H0 = W_MPO
        else
            H0 = +(H0, W_MPO; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
        end
    end
    return H0
end


function build_seed(sites, params)
    HS = nothing
    
    if !isnothing(params.S)
        println("Using seed for TCI: $(params.S)")
        _, S_MPO, _ = Quantics_TCI(params.S, Float64, sites, params.tci_tol)
        HS = S_MPO
    end
    return HS
end


# src/hamiltonians/mpo_construction.jl

function build_fock(sys::System)
    sites = sys.sites
    params = sys.params
    
    # 1. Get the bond-order coefficients via TCI
    fe = FockEvaluator1D(sys)
    _, F_MPO, _ = Quantics_TCI(x -> fe(x), Float64, sites, params.tci_tol)
    
    # 2. Re-use your hopping MPO logic
    # T_R and T_L are the bare sum_{i} c†_{i+1}c_i and sum_{i} c†_i c_{i+1}
    T_R, T_L = _build_translation_chain(sites) # Assuming this helper exists
    
    # 3. Apply the spatial modulation (The "Fock Hopping")
    # This creates sum_i F(i) * c†_{i+1}c_i
    VF_R = apply(F_MPO, T_R; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    
    # This creates sum_i F(i) * c†_i c_{i+1} (The h.c. part)
    # Since F(i) is real, dag(F_MPO) == F_MPO
    VF_L = apply(T_L, ITensors.dag(F_MPO); cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
    
    # Sum them to get the Hermitian Fock MPO
    return +(VF_R, VF_L; cutoff=params.itensors_tol, maxdim=params.itensors_maxdim)
end