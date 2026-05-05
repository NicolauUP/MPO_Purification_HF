
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