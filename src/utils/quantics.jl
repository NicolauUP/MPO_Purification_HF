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
        

        val_i = i - 1
        val_j = j - 1
        

        V = ITensor(1.0)
        
        for k in 1:L

            shift = L - k
            

            bit_i = (val_i >> shift) & 1
            bit_j = (val_j >> shift) & 1
            

            W_bra = bra_states[bit_i + 1][k]
            W_ket = ket_states[bit_j + 1][k]
            


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

function square_lattice_decoder(i::Integer, L::Integer)
    idx = i 
    
    x = zero(idx)
    y = zero(idx)
    
    for k in 0:(div(L,2)-1)
        # Even bits of idx (0, 2, 4...) control Y
        y |= ((idx >> (2k))     & 1) << k
        # Odd bits of idx (1, 3, 5...) control X
        x |= ((idx >> (2k + 1)) & 1) << k
    end
    
    return x, y
end
