

export System
struct System{Tt, Tu, Tw}
    d     :: Int
    n     :: Int
    sites :: Vector{Index{Int64}}
    t     :: Tt
    U     :: Tu
    W     :: Tw
end

function System(d::Int, n::Int, t::Tt, U::Tu, W::Tw, L::Int) where {Tt, Tu, Tw}
    if L % n != 0
        throw(ArgumentError("L=$L must be divisible by n=$n"))
    end
    sites = siteinds("Qubit", L)
    return System{Tt, Tu, Tw}(d, n, sites, t, U, W)
end


function Base.show(io::IO, sys::System)
    L = length(sys.sites)
    println(io, "System:")
    println(io, "  d     = $(sys.d)D")
    println(io, "  L     = $L sites (n=$(sys.n) sublattices)")
    println(io, "  t     :: $(typeof(sys.t))")
    println(io, "  U     :: $(typeof(sys.U))")
    println(io, "  W     :: $(typeof(sys.W))")
end