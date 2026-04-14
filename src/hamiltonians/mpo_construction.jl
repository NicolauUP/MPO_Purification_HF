# src/hamiltonian/mpo_construction.jl



function _build_translation_chain(sites)
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

function build_W(sites, params)
    if !isnothing(params.W)
        _, W_MPO, _ = Quantics_TCI(params.W, Float64, sites, params.tci_tol)
        return W_MPO
    else
        return nothing
    end
end


function build_H0(sites, params)
    T_R, T_L = _build_translation_chain(sites)
    H0 = nothing
    if params.t isa Number
        println("Using constant hopping t = $(params.t)")
        H0 = params.t * (T_R + T_L)
    elseif params.t isa Function
        println("Using function hopping t(x) = $(params.t)")
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


function build_seed(sites, params)
    HS = nothing
    
    if !isnothing(params.S)
        println("Using seed for TCI: $(params.S)")
        _, S_MPO, _ = Quantics_TCI(params.S, Float64, sites, params.tci_tol)
        HS = S_MPO
    end
    return HS
end