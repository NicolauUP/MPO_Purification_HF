
# -------------------------------------
# HELPERS 
# -------------------------------------
ITensors.op(::OpName"σ+", ::SiteType"Qubit") = [0 1; 0 0]
ITensors.op(::OpName"σ-", ::SiteType"Qubit") = [0 0; 1 0]
ITensors.op(::OpName"P+", ::SiteType"Qubit") = [0 0; 0 1]
ITensors.op(::OpName"P-", ::SiteType"Qubit") = [1 0; 0 0]
#ITensors.op(::OpName"Sz", ::SiteType"Qubit") = [1 0; 0 -1]

#= 
Remember to fix this eventualy! Sz or sigmaz can not be used as an operator name, since ITensors.jl already defines it for the spin-1/2 case.=# 

export Identity_MPO
function Identity_MPO(sites::Vector{<:Index})
    return MPO(sites, "Id")
end

"""
    CosineHopping(offset, amplitude, frequency, phase=0.0)

Callable nearest-neighbour hopping profile
`offset + amplitude * cos(frequency * x + phase)`. For a one-dimensional
QTT chain, this type additionally selects the exact finite-state MPO builder
instead of generic tensor-cross interpolation.
"""
struct CosineHopping <: Function
    offset::Float64
    amplitude::Float64
    frequency::Float64
    phase::Float64
    function CosineHopping(
        offset::Real,
        amplitude::Real,
        frequency::Real,
        phase::Real=0.0,
    )
        all(isfinite, (offset, amplitude, frequency, phase)) || throw(ArgumentError(
            "CosineHopping parameters must be finite",
        ))
        return new(
            Float64(offset), Float64(amplitude), Float64(frequency), Float64(phase),
        )
    end
end

(hopping::CosineHopping)(x::Real) =
    hopping.offset + hopping.amplitude * cos(hopping.frequency * Float64(x) + hopping.phase)
