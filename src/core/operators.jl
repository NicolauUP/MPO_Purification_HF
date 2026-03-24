
# -------------------------------------
# HELPERS 
# -------------------------------------
ITensors.op(::OpName"σ+", ::SiteType"Qubit") = [0 1; 0 0]
ITensors.op(::OpName"σ-", ::SiteType"Qubit") = [0 0; 1 0]
ITensors.op(::OpName"P+", ::SiteType"Qubit") = [0 0; 0 1]
ITensors.op(::OpName"P-", ::SiteType"Qubit") = [1 0; 0 0]
ITensors.op(::OpName"σz", ::SiteType"Qubit") = [1 0; 0 -1]

function Identity_MPO(sites::Vector{<:Index})
    L = length(sites)
    os = OpSum()
    for i in 1:L
        os += 1, "Id", i
    end
    return MPO(os, sites) / L
end