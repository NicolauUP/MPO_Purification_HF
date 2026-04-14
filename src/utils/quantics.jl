# src/utils/quantics.jl
function precompute_qtt_states(sites::Vector{<:Index})
    L = length(sites)
    

    bra_states = (
        [dag(state(sites[k]', "0")) for k in 1:L],
        [dag(state(sites[k]', "1")) for k in 1:L]
    )
    
    ket_states = (
        [state(sites[k], "0") for k in 1:L],
        [state(sites[k], "1") for k in 1:L]
    )
    
    return bra_states, ket_states
end

function MatrixChecker(mpo::MPO,sites::Vector{<:Index}, i::Int, j::Int, bra_states, ket_states)
    L = length(sites)
    
    # Convert 1-based grid index to 0-based binary value
    val_i = i - 1
    val_j = j - 1
    
    # Initialize the contraction accumulator as a scalar 1.0
    V = ITensor(1.0)
    
    for k in 1:L
        # MSB-first ordering: the shift amount decreases as we move down the chain
        shift = L - k
        
        # Hardware-level bit extraction (extracts a 0 or a 1)
        bit_i = (val_i >> shift) & 1
        bit_j = (val_j >> shift) & 1
        
        # Fetch precomputed tensors (+1 because Julia is 1-indexed)
        W_bra = bra_states[bit_i + 1][k]
        W_ket = ket_states[bit_j + 1][k]
        
        # Contract the physical indices to get the virtual matrix, 
        # then multiply it into the growing environment vector V.
        V *= (mpo[k] * W_bra * W_ket)
    end
    
    return scalar(V)
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