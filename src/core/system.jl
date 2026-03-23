

export System
"""
    System{T, V, W}

Holds all physical and numerical parameters for a spinless lattice HF calculation.

# Fields
- `L`     : total number of sites (= n * dL)
- `n`     : number of sublattice sites / orbitals per unit cell  
- `sites` : ITensors qubit site indices, consistent across all MPOs
- `t`     : hopping amplitude — scalar, Vector, or MPO (quasiperiodic)
- `U`     : interaction strength — scalar, Vector, or MPS (site-dependent)
- `mu`    : chemical potential — scalar, Vector, or MPS (site-dependent)
"""
struct System{T,V,W}
    d :: Int #spatial dimension (1D, 2D, etc.)
    n :: Int #number of orbitals/sublattices per site
    sites ::Vector{Index{Int64}}
    t :: T #hopping may be a scalar or a MPO!
    U :: V #interaction may be a scalar or a MPO!
    W :: W #external potential may be a scalar or a MPO!
end



"""
    System(L, n, t, U, mu)

Outer constructor — builds the ITensors qubit siteinds automatically.
`t`, `U`, `W` can each independently be a scalar, Vector, or MPS/MPO.
"""
function System(d::Int,n::Int, t::T, U::V, W::W, L::Int) where {T,V,W}
    if L % (n*d) != 0
        thow(ArgumentError("L=$L must be divisible by n=$n (number of sublattices) and d=$d (spatial dimension)"))
    end
    sites = siteinds("Qubit", L)
    return System{T,V,W}(L, n, sites, t, U, W)
end



"""
    num_cells(sys::System)

Returns the number of unit cells dL = L / n.
"""
num_cells(sys::System) = div(sys.L, sys.n)


"""
    Base.show(io::IO, sys::System)

Pretty-print a System in the REPL.
"""
function Base.show(io::IO, sys::System)
    L = sys.L
    println(io, "System:")
    println(io, " L    = $(sys.L) spin-sites")
    println(io, "t   = $(typeof(sys.t))")
    println(io, "U   = $(typeof(sys.U))")
    println(io, "W   = $(typeof(sys.W))")
end