# src/hamiltonian/mpo_construction.jl



function _build_translation_chain(sites; cutoff,maxdim)
    L = length(sites)

    T_R_opsum = OpSum()
    T_L_opsum = OpSum()

    for l in 1:L
        opsum_R_temp = OpSum()
        opsum_L_temp = OpSum()

        opsum_R_temp += 1, "σ+",l
        opsum_L_temp += 1, "σ-",l

        for m in l+1:L
            opsum_R_temp *= 1, "σ-",m
            opsum_L_temp *= 1, "σ+",m

        end

        T_R_opsum += opsum_R_temp
        T_L_opsum += opsum_L_temp

    end

    T_R_MPO = MPO(T_R_opsum, sites)
    T_L_MPO = MPO(T_L_opsum, sites)

return T_R_MPO, T_L_MPO
end

function build_H0(sites, params)
    T_R, T_L = _build_translation_chain(sites; cutoff=1e-10, maxdim=100)
    H0 = nothing
    if params.t isa Number
        H0 = params.t * (T_R + T_L)
    elseif params.t isa Function
        #= 
        TODO: Obtain the MPO of function t(x)
        =#
        T_MPO = Quantics_TCI(params.t, Float64, sites, 1e-6)[1] #[1] gets the MPO!
        H0 = apply(T_MPO, +(T_R, T_L; cutoff=cutoff, maxdim=maxdim), cutoff=cutoff, maxdim=maxdim)
    end

    if !isnothing(params.W)
        W_MPO = build_W(sites, params)
        H0 = +(H0, W_MPO; cutoff=cutoff, maxdim=maxdim)
    end
    return H0
end