module MPO_MeanField

using ITensors, ITensorMPS
using Random, SparseArrays
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf
using HDF5

include("core/operators.jl")
include("core/system.jl")
include("core/public_configuration.jl")
include("../src/hamiltonians/mpo_construction.jl")
include("../src/utils/quantics.jl")
include("../src/utils/runtime.jl")
include("../src/utils/profiling.jl")
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
include("../src/utils/qtt_charge_analysis.jl")
include("../src/io/checkpoint.jl")
include("kpm/public_kpm.jl")
include("kpm/qtt_mps_probe.jl")

export Parameters1D, ParametersSquare, ChainModel, SquareModel, GraphModel, CosineHopping, QTTSettings, SCFSettings, KPMSettings, RuntimeSettings, RuntimeReport, CaseSpec, CampaignSpec, PreflightIssue, PreflightReport, SolveResult, StoredResult, legacy_parameters, public_configuration, preflight, runtime_preflight, campaign_case, run_campaign, solve, write_result, read_result, write_mpo_checkpoint, read_mpo_checkpoint, qtt_bits, qtt_levels, System, PurificationResult, PurificationWorkStats, SCFIterationRecord, SCFDiagnostics, KPMIterationRecord, KPMDiagnostics, scf_diagnostics, extract_hartree_mpo_1d, extract_hartree_mpo_square, extract_hartree_mpo_tensorial_1d, extract_hartree_mpo_tensorial_square, extract_hartree_mpo_binary_carry_1d, extract_hartree_mpo_binary_carry_square, extract_hartree_mpo_binary_carry_square_adjacency, square_neighbour_adjacency_mpo, extract_fock_mpo_1d, extract_fock_mpo_binary_carry_1d, extract_fock_mpo_square_horizontal, extract_fock_mpo_square_vertical, extract_fock_mpo_binary_carry_square_horizontal, extract_fock_mpo_binary_carry_square_vertical, extract_mean_fields, run_scf!, MatrixChecker, construct_rho_0, perform_purification, perform_purification_sp2, build_translation_square, square_lattice_decoder, square_lattice_index, square_neighbours, square_undirected_bonds, validate_spectral_bounds, square_scf_spectral_bounds, verify_spectral_bounds_exact, nearest_neighbor_hf_energy_1d, nearest_neighbor_hf_energy_square, observables, observables_1d, observables_square, density_diagonal_mps, centered_density_mps, qtt_charge_ipr, qtt_multiscale_d2, qtt_fourier_square, qtt_fourier_amplitude, qtt_mps_amplitude, maybe_collect_garbage!
end # module MPO_MeanField
