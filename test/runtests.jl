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
include("test_basis_conventions.jl")
include("test_hamiltonian_1d.jl")
include("test_hamiltonian_square.jl")
include("test_fields_1d.jl")
include("test_energy_1d.jl")
include("test_square_scf_block.jl")
include("test_scf_noninteracting.jl")
include("test_scf_interacting_dense_reference.jl")
include("test_purification_result.jl")
include("test_spectral_bounds.jl")
include("test_purification_validation.jl")
include("test_sp2.jl")
include("test_sp2_polynomials.jl")
include("test_purification_crossvalidation.jl")
include("test_sp2_difficult_cases.jl")
include("test_mcweeny_mu.jl")
include("test_progress_output.jl")
