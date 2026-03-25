# src/utils/quantics.jl

"""
    Convert_To_Binary(n, NSites)

Convert an integer n to a binary string vector of length NSites.
Example: Convert_To_Binary(3, 4) = ["0", "0", "1", "1"]
"""
function Convert_To_Binary(n::Int64, NSites::Int64)
    bits = last(bitstring(n), NSites)
    VectorString = String[]
    for i in bits
        push!(VectorString, string(i))
    end
    return VectorString
end

"""
    BasisStateMPS(n, sites)

Create the computational basis MPS |n⟩ for integer n
in the binary representation of the Hilbert space.
"""
function BasisStateMPS(n::Int64, sites::Vector{<:Index})
    return MPS(sites, Convert_To_Binary(n, length(sites)))
end

"""
    MatrixChecker(mpo, sites, i, j)

Evaluate the matrix element ⟨i|MPO|j⟩ between computational
basis states |i⟩ and |j⟩ specified by their integer indices.
Useful for debugging and verifying MPO constructions.
"""
function MatrixChecker(mpo::MPO, sites::Vector{<:Index}, i::Int64, j::Int64)
    Ψ_i = BasisStateMPS(i, sites)
    Ψ_j = BasisStateMPS(j, sites)
    return inner(Ψ_i', mpo, Ψ_j)
end

"""
    Quantics_TCI(f, eltype, sites, ϵ)

Compute the quantics TCI of a function f on the integer domain {0, ..., 2^L - 1}.
Returns (QTT, mpo, mps) where mpo is the diagonal MPO representation of f.
"""
function Quantics_TCI(f::Function, eltype::Type{<:Number}, sites::Vector{<:Index}, ϵ::Float64)
    NSites = length(sites)
    XVals = range(0, 2^NSites - 1; length=2^NSites)
    QTT, _, _ = QuanticsTCI.quanticscrossinterpolate(eltype, f, XVals; tolerance=ϵ)
    TT = TCI.tensortrain(QTT.tci)
    mps = MPS(TT; sites)
    mpo = MPO(NSites)
    for i in 1:NSites
        mpo.data[i] = Quantics._asdiagonal(mps.data[i], sites[i])
    end
    return QTT, mpo, mps
end

"""
    Quantics_TCI(f, eltype, sites, xmin, xmax, ϵ)

Compute the quantics TCI of a function f on the physical domain [xmin, xmax]
mapped onto 2^L points. Returns (QTT, mpo, mps).
"""
function Quantics_TCI(f::Function, eltype::Type{<:Number},
                      sites::Vector{<:Index}, xmin::Float64, xmax::Float64, 
                      ϵ::Float64)
    NSites = length(sites)
    XVals = range(xmin, xmax; length=2^NSites)
    QTT, _, _ = QuanticsTCI.quanticscrossinterpolate(eltype, f, XVals; tolerance=ϵ)
    TT = TCI.tensortrain(QTT.tci)
    mps = MPS(TT; sites)
    mpo = MPO(NSites)
    for i in 1:NSites
        mpo.data[i] = Quantics._asdiagonal(mps.data[i], sites[i])
    end
    return QTT, mpo, mps
end