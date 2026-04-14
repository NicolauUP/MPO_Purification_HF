

struct ModelParameters{Tt, Tu, Tw, Ts}
    L::Int
    t::Tt #Type of the hopping, maybe either a number or a function to be used with TCI
    U::Tu #Type of the interaction
    W::Tw #Type of the potential.
    S::Ts #Type of the seed for TCI
    tci_tol::Float64
    itensors_tol::Float64
    itensors_maxdim::Int
    density::Float64
    purification_steps::Int
    scf_mixing::Float64
    scf_tol::Float64
    scf_max_iterations::Int
end

mutable struct System{P}
    params::P
    sites::Vector{Index{Int64}}
    H0::MPO # This is just H0 + W
    VH::MPO # Dynamic: Hartree Potential
    VF::MPO # Dynamic: Fock Potential
    ρ::MPO # Dynamic: Density Matrix

    bra_states::Any
    ket_states::Any
end


function System(params::ModelParameters)
    sites = ITensors.siteinds("Qubit", params.L)
    

    H_static = build_H0(sites, params)

    VH_init = build_seed(sites, params)
    VF_init = Identity_MPO(sites) * 0.0 # We start with no Fock potential, but we could also build a seed for it.
    rho_init = Identity_MPO(sites) * 0.0 #nothing

    bra, ket = precompute_qtt_states(sites)
    
    return System{typeof(params)}(
        params,
        sites,
        H_static,
        VH_init,
        VF_init,
        rho_init,
        bra,
        ket
    )
end

function Base.show(io::IO, sys::System)
    println(io, "System (L=$(sys.params.L))")
    println(io, "  ├─ Hopping (t): $(typeof(sys.params.t))")
    println(io, "  ├─ Interaction (U): $(typeof(sys.params.U))")
    println(io, "  ├─ Potential (W): $(isnothing(sys.params.W) ? "None" : typeof(sys.params.W))")
    println(io, "  ├─ TCI Precision: $(sys.params.tci_tol)")
    println(io, "  ├─ ITensors Precision: $(sys.params.itensors_tol)")
    println(io, "  └─ ITensors MaxDim: $(sys.params.itensors_maxdim)")
    println(io, "  [Dynamic State]")
    println(io, "  ├─ VH MaxLinkDim: $(maxlinkdim(sys.VH))")
    println(io, "  └─ ρ MaxLinkDim:  $(maxlinkdim(sys.ρ))")
end
