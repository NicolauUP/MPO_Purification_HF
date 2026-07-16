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

    """
        MatrixChecker(mpo, sites, i, j, bra_states, ket_states)

    Return the one-based matrix element `⟨i-1|mpo|j-1⟩`. Binary basis
    labels are decoded most-significant bit first: `sites[1]` carries the
    highest-order bit and `sites[end]` the lowest-order bit.
    """
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

function zero_mpo(sites::Vector{<:Index})
    return 0.0 * Identity_MPO(sites)
end

function diagonal_mpo_from_function(
    f::Function,
    eltype::Type{<:Number},
    sites::Vector{<:Index},
    tolerance::Float64,
)
    if length(sites) == 1
        diagonal = OpSum()
        diagonal += convert(eltype, f(0)), "P-", 1
        diagonal += convert(eltype, f(1)), "P+", 1
        return MPO(diagonal, sites)
    end

    try
        _, mpo, _ = Quantics_TCI(f, eltype, sites, tolerance)
        return mpo
    catch err
        if err isa ErrorException && occursin("maxsamplevalue is zero", err.msg)
            return zero_mpo(sites)
        end
        rethrow()
    end
end

function safe_relative_change(diff_norm::Real, reference_norm::Real; floor::Real=sqrt(eps(Float64)))
    diff_norm >= 0 || throw(ArgumentError("diff_norm must be nonnegative"))
    reference_norm >= 0 || throw(ArgumentError("reference_norm must be nonnegative"))
    floor > 0 || throw(ArgumentError("floor must be positive"))
    return diff_norm / max(reference_norm, floor)
end

"""
    print_iteration_progress(io, label, iteration, total, details; overwrite=io isa Base.TTY)

Print a single iteration status. On an interactive terminal, the next status
replaces this one; redirected output remains newline-delimited for logs.
"""
function print_iteration_progress(
    io::IO,
    label::AbstractString,
    iteration::Integer,
    total::Integer,
    details::AbstractString;
    overwrite::Bool=io isa Base.TTY,
)
    message = "$label $iteration/$total | $details"
    if overwrite
        print(io, "\r\e[2K", message)
        flush(io)
    else
        println(io, message)
    end
    return overwrite
end

function finish_iteration_progress(io::IO, overwrite::Bool)
    overwrite && println(io)
    return nothing
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

function square_lattice_decoder(i::Number, L::Integer)
    idx = Int(i)
    
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
