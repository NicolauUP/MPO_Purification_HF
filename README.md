# HartreeFockMPO.jl

A Julia package for Hartree-Fock mean-field theory of short-range interacting
spinless lattice systems, using Matrix Product Operator (MPO) techniques.

## Physical Setup

- Spinless fermions on a lattice with `2^L` sites in binary (qubit) representation
- Short-range (first-neighbor) hopping, Hubbard-like on-site interactions
- Support for uniform or spatially modulated (quasiperiodic) hopping and static potential
- Modulations specified as plain Julia functions — TCI conversion handled internally
- Initial state symmetry breaking via TCI-generated seed potentials (e.g., to induce Charge Density Waves)
- Canonical Hartree-Fock via MPO purification and Tensor Cross Interpolation (TCI)

## Method Overview

1. **Static Hamiltonian & Seed Construction**: The static single-particle Hamiltonian (`H0`) is written as an MPO using a binary representation of the lattice. If `t` or `W` are functions, they are converted to diagonal MPOs via Quantics TCI. An initial symmetry-breaking potential (`seed_potential`) is also constructed to jumpstart the SCF loop.

2. **McWeeny Purification**: The effective Hamiltonian (`H_eff = H0 + VH`) is purified into the ground state density matrix `ρ` using an adaptive scheme. It starts with a trace-correcting linear update and switches to McWeeny (3ρ² - 2ρ³) when idempotency is close.

3. **TCI Extraction**: Tensor Cross Interpolation (TCI) is used dynamically to:
   - Extract the Hartree potential (diagonal MPO) directly from the current density matrix by sampling `⟨i|ρ|i⟩`.
   - *[Upcoming]* Extract the Fock exchange terms (off-diagonal).

4. **SCF Loop**: The newly extracted Hartree MPO is mixed with the previous iteration's potential (Density Mixing) and fed back into the effective Hamiltonian. The system is iteratively re-purified until self-consistency is reached.

## Package Structure

```text
HartreeFockMPO/
├── Project.toml                    # Package dependencies
├── README.md                       # This file
├── test_main.jl                    # Main test script
│
└── src/
    ├── HartreeFockMPO.jl           # Module entry point — includes + exports
    │
    ├── core/
    │   ├── operators.jl            # Custom ITensors qubit operator definitions
    │   │                           # (σ+, σ-, P+, P-, σz) + Identity_MPO
    │   └── system.jl               # ModelParameters and mutable System structs
    │                               # Initializes static H0 and dynamic state (VH, ρ)
    │
    ├── hamiltonians/
    │   └── mpo_construction.jl     # Build H0 and Seeds as MPOs
    │                               # _build_translation_chain: bare hopping MPOs
    │                               # build_H0: Kinetic + static potential W
    │                               # build_seed: TCI conversion of seed function
    │
    ├── purification/
    │   └── mcweeny.jl              # Density matrix via purification
    │                               # construct_rho_0(sys, params, H_max, H_min)
    │                               #   - Bakes H_eff = H0 + VH into initial guess
    │                               # perform_purification(ρ0, params; verbose)
    │                               #   - Trace-correcting -> McWeeny (3ρ²-2ρ³)
    │
    ├── tci/
    │   ├── modulations.jl          # Coordinate mappers & static modulations
    │   │                           #     index_to_xy_snake, etc.
    │   └── density_matrix.jl       # TCI extraction of Hartree/Fock from ρ
    │                               # build_hartree(sys): Extracts V_H MPO dynamically
    │                               # HartreeEvaluator1D: Maps TCI queries to ⟨j|ρ|j⟩
    │
    ├── hf/
    │   └── selfconsistent.jl       # SCF loop orchestrator
    │                               # run_scf!(sys, ...): Drives the purification ->
    │                               # extraction -> mixing cycle until convergence
    │
    └── utils/
        ├── quantics.jl             # TCI/Quantics utilities
        │                           # MatrixChecker: evaluates ⟨i|MPO|j⟩
        │                           # precompute_qtt_states: speeds up MatrixChecker
        └── observables.jl          # Physical observables
                                    # [TO BE IMPLEMENTED]

## Data Flow

```text
ModelParameters{Tt, Tu, Tw, Ts}
  │
  └── System(params)
        ├── sites, bra_states, ket_states
        ├── H0 = build_H0(sites, params)            # Static environment
        ├── VH = build_seed(sites, params)          # Symmetry-breaking seed
        └── rho = Id                                # Dummy init

[SCF LOOP: run_scf!(sys, H_max, H_min)]
  │
  ├──► ρ0 = construct_rho_0(sys, ...)          # Bakes H_eff = H0 + VH into initial guess
  │
  ├──► ρ = perform_purification(ρ0, params)    # Pushes eigenvalues to {0, 1}
  │
  ├──► new_VH = build_hartree(sys)             # TCI probes new ρ, builds next Hartree MPO
  │
  └──► sys.VH = α * new_VH + (1 - α) * sys.VH  # Mix potentials to stabilize CDW / loop

  ## Key Structs

```julia
ModelParameters{Tt, Tu, Tw, Ts}
  L                 :: Int       # number of qubits (lattice has 2^L sites)
  t                 :: Tt        # hopping:     Number | Function
  U                 :: Tu        # interaction: Number | Function
  W                 :: Tw        # static potential: Nothing | Function
  seed_potential    :: Ts        # initial symmetry breaking: Nothing | Function
  tci_tol           :: Float64   # TCI tolerance for function → MPO
  itensors_tol      :: Float64   # ITensors SVD cutoff
  itensors_maxdim   :: Int       # ITensors max bond dimension
  density           :: Float64   # filling fraction (Ne = 2^L * density)
  purification_steps:: Int       # max McWeeny iterations

mutable struct System{P}
  params     :: P                    # ModelParameters
  sites      :: Vector{Index{Int64}} # ITensors qubit indices, length L
  H0         :: MPO                  # Static: Kinetic + W
  VH         :: MPO                  # Dynamic: Seed -> Hartree Potential
  rho        :: MPO                  # Dynamic: Density Matrix
  bra_states :: Any                  # Precomputed bases for MatrixChecker
  ket_states :: Any                  # Precomputed bases for MatrixChecker