using Test
using Random
using MPO_MeanField
using ITensors, ITensorMPS
using LinearAlgebra

Random.seed!(0x5eed)

function dense_matrix(mpo, sys)
    N = 2^sys.params.L
    return [MatrixChecker(mpo, sys.sites, i, j, sys.bra_states, sys.ket_states) for i in 1:N, j in 1:N]
end

include("test_construction.jl")
include("test_zero_safe.jl")
