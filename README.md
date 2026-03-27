# HartreeFockMPO.jl

A Julia package for Hartree-Fock mean-field theory of short-range interacting
spinless lattice systems, using Matrix Product Operator (MPO) techniques.

## Physical Setup

- Spinless fermions on a lattice with `2^L` sites in binary (qubit) representation
- Short-range (first-neighbor) hopping, Hubbard-like on-site interactions
- Support for uniform or spatially modulated (quasiperiodic) hopping and potential
- Modulations specified as plain Julia functions — TCI conversion handled internally
- Canonical Hartree-Fock via MPO purification and Tensor Cross Interpolation (TCI)

## Method Overview

1. **H0 construction**: The single-particle Hamiltonian is written as an MPO
   using a binary (qubit) representation of the lattice. Hopping terms correspond
   to binary addition operators (σ+/σ- chains). If `t` or `W` are Julia functions,
   they are automatically converted to diagonal MPOs via Quantics TCI internally.

2. **McWeeny Purification**: The Hamiltonian MPO is purified to the ground state
   density matrix ρ using an adaptive scheme that starts with a trace-correcting
   linear update and switches to McWeeny (3P² - 2P³) when idempotency is close.

3. **TCI Extraction**: Tensor Cross Interpolation (TCI) is used to:
   - Build diagonal MPOs for spatially modulated hoppings and potentials
   - Extract Hartree (diagonal) and Fock (off-diagonal) terms from ρ

4. **SCF Loop**: The Hartree and Fock MPOs are fed back into H0 to form H_HF,
   which is re-purified until self-consistency is reached.

## Package Structure

```
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
    │   └── system.jl               # ModelParameters and System structs
    │                               # ModelParameters{Tt, Tu, Tw}:
    │                               #   L, t, U, W, tci_tol, itensors_tol,
    │                               #   itensors_maxdim, density, purification_steps
    │                               # System{P}:
    │                               #   params, sites, H0
    │                               #   H0 is built automatically via System(params)
    │
    ├── hamiltonians/
    │   └── mpo_construction.jl     # Build H0 as MPO
    │                               # _build_translation_chain: bare hopping MPOs
    │                               # build_W: TCI conversion of W::Function → MPO
    │                               # build_H0: dispatches on t type
    │                               #   - t::Number   → uniform hopping
    │                               #   - t::Function → TCI modulated hopping
    │                               # W added automatically if !isnothing(params.W)
    │
    ├── purification/
    │   └── mcweeny.jl              # Density matrix via purification
    │                               # construct_rho_0(sys, params, H_max, H_min)
    │                               #   N  = 2^L (Hilbert space = lattice sites)
    │                               #   Ne = round(Int, N * density)
    │                               #   Id built internally from sys.sites
    │                               # perform_purification(ρ0, params; verbose)
    │                               #   - trace-correcting update until cn ≈ 0.5
    │                               #   - switches to McWeeny (3P²-2P³)
    │                               #   - @warn on stuck or non-converged
    │
    ├── tci/
    │   ├── modulations.jl          # Model-specific modulation MPOs via TCI
    │   │                           # 1D: build_pi_modulation
    │   │                           #     build_quasiperiodic_modulation
    │   │                           # 2D: build_pi_modulation_2D
    │   │                           #     build_quasiperiodic_modulation_2D
    │   │                           # Coordinate mappers:
    │   │                           #     index_to_xy_lexi
    │   │                           #     index_to_xy_snake
    │   └── density_matrix.jl       # TCI extraction of Hartree/Fock from ρ
    │                               # [TO BE IMPLEMENTED]
    │
    ├── hf/
    │   └── selfconsistent.jl       # SCF loop
    │                               # H_HF = H0 + V_hartree + V_fock
    │                               # iterates until convergence
    │                               # [TO BE IMPLEMENTED]
    │
    └── utils/
        ├── quantics.jl             # TCI/Quantics utilities
        │                           # Convert_To_Binary: Int → binary string vector
        │                           # BasisStateMPS: Int → computational basis |n⟩
        │                           # MatrixChecker: evaluates ⟨i|MPO|j⟩
        │                           # Quantics_TCI: f::Function → diagonal MPO
        │                           #   - integer domain {0,...,2^L-1}
        │                           #   - physical domain [xmin, xmax]
        └── observables.jl          # Physical observables
                                    # [TO BE IMPLEMENTED]
```

## Data Flow

```
ModelParameters{Tt, Tu, Tw}
  │
  └── System(params)
        ├── sites = siteinds("Qubit", L)
        └── H0 = build_H0(sites, params)
                  ├── [t::Number]   → t * (T_R + T_L)
                  ├── [t::Function] → TCI(t) ⊗ (T_R + T_L)
                  └── [W::Function] → H0 + TCI(W)

ρ0 = construct_rho_0(sys, params, H_min, H_max)
      ├── N  = 2^L  (number of lattice sites)
      ├── Ne = round(Int, N * density)
      └── ρ0 = coeff_I * Id + coeff_H * H0

ρ = perform_purification(ρ0, params; verbose)
      ├── adaptive trace-correcting update
      └── switches to McWeeny when cn ≈ 0.5

[TO BE IMPLEMENTED]
V_hartree, V_fock = extract_hf_terms(ρ, sys)   # tci/density_matrix
H_HF = H0 + V_hartree + V_fock                 # hf/selfconsistent
loop until convergence
```

## Key Structs

```julia
ModelParameters{Tt, Tu, Tw}
  L                 :: Int       # number of qubits (lattice has 2^L sites)
  t                 :: Tt        # hopping:     Number | Function
  U                 :: Tu        # interaction: Number | Function
  W                 :: Tw        # potential:   Nothing | Function
  tci_tol           :: Float64   # TCI tolerance for function → MPO
  itensors_tol      :: Float64   # ITensors SVD cutoff
  itensors_maxdim   :: Int       # ITensors max bond dimension
  density           :: Float64   # filling fraction (Ne = 2^L * density)
  purification_steps:: Int       # max McWeeny iterations

System{P}
  params  :: P                       # ModelParameters
  sites   :: Vector{Index{Int64}}    # ITensors qubit indices, length L
  H0      :: MPO                     # single-particle Hamiltonian (+ potential)
```

## Dependencies

```julia
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using TCIITensorConversion
using LinearAlgebra
```

## Usage Example

```julia
# 1D chain, L=6 qubits (2^6 = 64 lattice sites), half filling
L    = 6
t    = -1.0
U    = 0.0
W(x) = 0.5 * cos(π * x)   # π-modulated potential, plain Julia function

params = ModelParameters(L, t, U, W, 1e-6, 1e-10, 100, 0.5, 40)
sys    = System(params)    # builds sites and H0 automatically

# Purify to density matrix
ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
ρ  = perform_purification(ρ0, params; verbose=1)

# Quasiperiodic hopping (golden ratio)
β       = (1 + √5) / 2
t_qp(x) = cos(2π * β * x)
params_qp = ModelParameters(L, t_qp, U, nothing, 1e-6, 1e-10, 100, 0.5, 40)
sys_qp    = System(params_qp)
```

## Status

| Component                           | Status             |
|-------------------------------------|--------------------|
| `core/operators.jl`                 | ✅ Done            |
| `core/system.jl`                    | ✅ Done            |
| `hamiltonians/mpo_construction.jl`  | ✅ Done (1D chain) |
| `purification/mcweeny.jl`           | ✅ Done            |
| `utils/quantics.jl`                 | ✅ Done            |
| `tci/modulations.jl`                | ✅ Done (1D)  |
| `tci/density_matrix.jl`             | 🔲 To implement    |
| `hf/selfconsistent.jl`              | 🔲 To implement    |
| `utils/observables.jl`              | 🔲 To implement    |
