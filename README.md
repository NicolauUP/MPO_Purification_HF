# HartreeFockMPO.jl

A Julia package for Hartree-Fock mean-field theory of short-range interacting 
spinless lattice systems, using Matrix Product Operator (MPO) techniques.

## Physical Setup

- Spinless fermions on a lattice with short-range (first-neighbor) hopping
- Hubbard-like on-site interactions
- Support for uniform, staggered, and quasiperiodic modulations of hopping and potential
- Canonical Hartree-Fock via MPO purification and Tensor Cross Interpolation (TCI)

## Method Overview

1. **H0 construction**: The single-particle Hamiltonian is written as an MPO
   using a binary (qubit) representation of the lattice. Hopping terms correspond
   to binary addition operators (σ+/σ- chains). Geometry is hardcoded per lattice type.

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
    │   └── system.jl               # System{Tt,Tu,Tw} struct
    │                               # Holds: d, n, sites, t, U, W
    │                               # t/U/W can be scalar, Vector, or MPO
    │
    ├── hamiltonian/
    │   └── mpo_construction.jl     # Build H0 as MPO
    │                               # _build_translation_chain: bare hopping MPO
    │                               # build_hopping_chain: dispatches on t type
    │                               #   - t::Float64 → uniform hopping
    │                               #   - t::MPO     → modulated hopping
    │                               # build_H0_chain: adds W potential
    │                               #   - W::Nothing → hopping only
    │                               #   - W::MPO     → hopping + potential
    │
    ├── purification/
    │   └── mcweeny.jl              # Density matrix via purification
    │                               # construct_rho_0: initial guess ρ0
    │                               #   inputs: sys, H, ϵ, maxχ, H_max, H_min, Ne
    │                               #   H_max, H_min, Ne are user-provided
    │                               # perform_purification: adaptive McWeeny loop
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
        │                           # Quantics_TCI: f → diagonal MPO (two variants)
        │                           #   - integer domain {0,...,2^L-1}
        │                           #   - physical domain [xmin, xmax]
        └── observables.jl          # Physical observables
                                    # [TO BE IMPLEMENTED]
```

## Data Flow

```
System{Tt, Tu, Tw}
  │
  ├── [if Tw=MPO] tci/modulations → W_mpo
  ├── [if Tt=MPO] tci/modulations → t_mpo
  │
  ├── hamiltonian/mpo_construction → H0
  │
  └── hf/selfconsistent:
        ρ0 = construct_rho_0(sys, H_HF, ϵ, maxχ, H_max, H_min, Ne)
        loop:
          ρ  = perform_purification(ρ0; maxχ, ϵ, max_steps, verbose)
          V_hartree, V_fock = extract_hf_terms(ρ, sys)   # tci/density_matrix
          H_HF = H0 + V_hartree + V_fock
        until convergence
```

## System Struct

```julia
System{Tt, Tu, Tw}
  d     :: Int                    # spatial dimensionality (1, 2, ...)
  n     :: Int                    # number of sublattices/orbitals
  sites :: Vector{Index{Int64}}   # ITensors qubit indices, length = L (total sites)
  t     :: Tt                     # hopping:     Float64 | MPO
  U     :: Tu                     # interaction: Float64 | MPO
  W     :: Tw                     # potential:   Nothing | MPO
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
using HartreeFockMPO

# Define a 1D chain, L=10, uniform hopping, no potential
sys = System(1, 1, 1.0, 0.0, nothing, 10)

# Build H0
H0 = build_H0_chain(sys; cutoff=1e-10, maxdim=100)

# Purify to density matrix at half filling
ρ0 = construct_rho_0(sys, H0, 1e-10, 100, 2.0, -2.0, 5)
ρ  = perform_purification(ρ0; maxχ=100, ϵ=1e-10, verbose=1)

# Build a quasiperiodic potential (golden ratio)
β = (1 + √5) / 2
W = build_quasiperiodic_modulation(sys, 0.5, β, 0.0)

# Rebuild system with potential and run SCF
sys_W = System(1, 1, 1.0, 1.0, W, 10)
H0_W  = build_H0_chain(sys_W; cutoff=1e-10, maxdim=100)
```

## Status

| Component                        | Status            |
|----------------------------------|-------------------|
| `core/operators.jl`              | ✅ Done           |
| `core/system.jl`                 | ✅ Done           |
| `hamiltonian/mpo_construction.jl`| ✅ Done (1D chain)|
| `purification/mcweeny.jl`        | ✅ Done           |
| `utils/quantics.jl`              | ✅ Done           |
| `tci/modulations.jl`             | ✅ Done (1D + 2D) |
| `tci/density_matrix.jl`          | 🔲 To implement   |
| `hf/selfconsistent.jl`           | 🔲 To implement   |
| `utils/observables.jl`           | 🔲 To implement   |
