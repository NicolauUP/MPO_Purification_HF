

struct ModelParameters{Tt, Tu, Tw}
    L::Int
    t::Tt #Type of the hopping, maybe either a number or a function to be used with TCI
    U::Tu #Type of the interaction
    W::Tw #Type of the potential.
end

struct System{P}
    params::P
    sites::Vector{Index{Int64}}
    H0::MPO
    W::Union{Nothing, MPO}
    H_static::MPO # This is just H0 + W
end


function System(params)
    sites = ITensors.siteinds("Qubit", params.L)
    H0 = MPO(sites,"Id")
    W = isnothing(params.W)  ? nothing : build_W(sites, params)
    H_static = isnothing(params.W) ? H0 : H0 + W  
    return System(params, sites, H0, W, H_static)
end





function Base.show(io::IO, sys::System)
    println(io, "System:")
    println(io, "L: ", sys.params.L)
    println(io, "Type of t: ", typeof(sys.params.t))
    println(io, "Type of U: ", typeof(sys.params.U))
    println(io, "Type of W: ", typeof(sys.params.W))
end

#=

TODO: 
    - Add a method for Base.show to print the system
=#