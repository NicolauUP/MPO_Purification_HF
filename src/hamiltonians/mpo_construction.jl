# src/hamiltonian/mpo_construction.jl

export build_hopping_chain
# -------------------------------------
# 1D Chain - uniform hopping t::Float64
# -------------------------------------

function build_hopping_chain(sys::System{Float64, V, W}; cutoff=1e-10, maxdim=100) where {V,W}
    sites = sys.sites
    L = length(sites)
    t = sys.t

    T_R_opsum = OpSum()
    T_L_opsum = OpSum()

    for l in 1:L
        opsum_R_temp = OpSum()
        opsum_L_temp = OpSum()

        opsum_R_temp += t, "σ+",l
        opsum_R_temp += t, "σ-",L

        for m in l+1:L
            opsum_R_temp *= t, "σ-",m
            opsum_L_temp *= t, "σ+",m

        end

        T_R_opsum += opsum_R_temp
        T_L_opsum += opsum_L_temp

    end

    T_R_MPO = MPO(T_R_opsum, sites)
    T_L_MPO = MPO(T_L_opsum, sites)

return +(T_R_MPO, T_L_MPO; cutoff=cutoff, maxdim=maxdim)
end


# -------------------------------------
# 1D Chain - modulated hopping t::MPO 
# -------------------------------------

function build_hopping_chain(sys::System{MPO, V, W}; cutoff=1e-10, maxdim=100) where {V,W}
    sites = sys.sites
    L = length(sites)
    t_mpo = sys.t

    T_R_opsum = OpSum()
    T_L_opsum = OpSum()

    for l in 1:L
        opsum_R_temp = OpSum()
        opsum_L_temp = OpSum()

        opsum_R_temp += 1, "σ+",l
        opsum_R_temp += 1, "σ-",L

        for m in l+1:L
            opsum_R_temp *= 1, "σ-",m
            opsum_L_temp *= 1, "σ+",m

        end

        T_R_opsum += opsum_R_temp
        T_L_opsum += opsum_L_temp

    end

    T_R_MPO = MPO(T_R_opsum, sites)
    T_L_MPO = MPO(T_L_opsum, sites)

    T_R_mod = apply(t_mpo, T_R_MPO; cutoff=cutoff, maxdim=maxdim)
    T_L_mod = apply(t_mpo, T_L_MPO; cutoff=cutoff, maxdim=maxdim)
return +(T_R_MPO, T_L_MPO; cutoff=cutoff, maxdim=maxdim)
end