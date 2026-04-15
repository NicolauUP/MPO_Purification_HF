module MPO_MeanField

using ITensors, ITensorMPS
using Quantics, QuanticsTCI
using TCIITensorConversion
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf

include("../src/core/operators.jl")
include("../src/core/system.jl")
include("../src/hamiltonians/mpo_construction.jl")
include("../src/utils/quantics.jl")
include("../src/purification/mcweeny.jl")
include("../src/tci/modulations.jl")
include("../src/tci/density_matrix.jl")
include("../src/hf/self_consistent.jl")

export ModelParameters, System, extract_hartree_mpo_1d, extract_fock_mpo_1d, run_scf!, MatrixChecker, construct_rho_0, perform_purification
end # module MPO_MeanField
