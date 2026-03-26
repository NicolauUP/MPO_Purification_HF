

struct ModelParameters{Tt, Tu, Tw}
    L::Int
    t::Tt #Type of the hopping, maybe either a number or a function to be used with TCI
    U::Tu #Type of the interaction
    W::Tw #Type of the potential.
    tci_tol::Float64
    itensors_tol::Float64
    itensors_maxdim::Int
    density::Float64
    purification_steps::Int
end

struct System{P}
    params::P
    sites::Vector{Index{Int64}}
    H0::MPO # This is just H0 + W
end

function System(params::ModelParameters)
    sites = ITensors.siteinds("Qubit", params.L)
    

    H_static = build_H0(sites, params)
    
    return System{typeof(params)}(params, sites, H_static)
end


function Base.show(io::IO, sys::System)
    println(io, "System (L=$(sys.params.L))")
    println(io, "  ├─ Hopping (t): $(typeof(sys.params.t))")
    println(io, "  ├─ Interaction (U): $(typeof(sys.params.U))")
    println(io, "  ├─ Potential (W): $(isnothing(sys.params.W) ? "None" : typeof(sys.params.W))")
    println(io, "  └─ TCI Precision: $(sys.params.tci_tol)")
    println(io, "  ├─ ITensors Precision: $(sys.params.itensors_tol)")
    println(io, "  └─ ITensors MaxDim: $(sys.params.itensors_maxdim)")
end
