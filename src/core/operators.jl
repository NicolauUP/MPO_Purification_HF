
# -------------------------------------
# HELPERS 
# -------------------------------------
ITensors.op(::OpName"σ+", ::SiteType"Qubit") = [0 1; 0 0]
ITensors.op(::OpName"σ-", ::SiteType"Qubit") = [0 0; 1 0]
ITensors.op(::OpName"P+", ::SiteType"Qubit") = [0 0; 0 1]
ITensors.op(::OpName"P-", ::SiteType"Qubit") = [1 0; 0 0]
ITensors.op(::OpName"σz", ::SiteType"Qubit") = [1 0; 0 -1]

export Identity_MPO
function Identity_MPO(sites::Vector{<:Index})
    return MPO(sites, "Id")
end