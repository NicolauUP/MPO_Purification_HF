module MPO_MeanField

using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf

include("core/operators.jl")
include("core/system.jl")
include("../src/hamiltonians/mpo_construction.jl")
include("../src/utils/quantics.jl")
include("../src/utils/runtime.jl")
include("../src/purification/result.jl")
include("../src/purification/initialization.jl")
include("../src/purification/sp2.jl")
include("../src/purification/mcweeny.jl")
include("../src/purification/palser_manolopoulos.jl")
include("../src/purification/purification.jl")
include("../src/tci/modulations.jl")
include("../src/tci/density_matrix.jl")
include("../src/hf/self_consistent.jl")
include("../src/utils/observables.jl")

export Parameters1D, ParametersSquare, System, PurificationResult, PurificationWorkStats, extract_hartree_mpo_1d, extract_hartree_mpo_square, extract_fock_mpo_1d, extract_fock_mpo_square_horizontal, extract_fock_mpo_square_vertical, extract_mean_fields, run_scf!, MatrixChecker, construct_rho_0, perform_purification, perform_purification_sp2, build_translation_square, square_lattice_decoder, square_lattice_index, square_neighbours, square_undirected_bonds, validate_spectral_bounds, verify_spectral_bounds_exact, nearest_neighbor_hf_energy_1d, nearest_neighbor_hf_energy_square, observables_1d, maybe_collect_garbage!
end # module MPO_MeanField
