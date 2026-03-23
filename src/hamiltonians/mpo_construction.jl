# src/hamiltonian/mpo_construction.jl

export build_hopping_chain

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
# -------------------------------------
# 1D Chain - uniform hopping t::Float64
# -------------------------------------

function build_hopping_chain(sys::System{Float64, V, W}; cutoff=1e-10, maxdim=100) where {V,W}
   T_R, T_L = _build_translation_chain(sys.sites; cutoff=cutoff, maxdim=maxdim)
   return +(sys.t * T_R, sys.t * T_L; cutoff=cutoff, maxdim=maxdim)
end


# -------------------------------------
# 1D Chain - modulated hopping t::MPO 
# -------------------------------------

function build_hopping_chain(sys::System{MPO, V, W}; cutoff=1e-10, maxdim=100) where {V,W}
    T_R, T_L = _build_translation_chain(sys.sites; cutoff=cutoff, maxdim=maxdim)
    T_R_mod = apply(sys.t, T_R; cutoff=cutoff, maxdim=maxdim)
    T_L_mod = apply(sys.t, T_L; cutoff=cutoff, maxdim=maxdim)
    return +(T_R_mod, T_L_mod; cutoff=cutoff, maxdim=maxdim)
end



function build_H0_chain(sys::System{T,V,Nothing}; cutoff=1e-10, maxdim=100) where {T,V}
    return build_hopping_chain(sys; cutoff=cutoff, maxdim=maxdim)
end

function build_H0_chain(sys::System{T,V,MPO}; cutoff=1e-10, maxdim=100) where {T,V}
    H_hop = build_hopping_chain(sys; cutoff=cutoff, maxdim=maxdim)
    return +(H_hop,  sys.W; cutoff=cutoff, maxdim=maxdim)
end