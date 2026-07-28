# MPO_MeanField.jl

`MPO_MeanField.jl` is a Julia research package for zero-temperature,
spinless, nearest-neighbour Hartree–Fock calculations on binary-encoded
lattices. It represents one-body Hamiltonians, density matrices, and
mean-field potentials as matrix-product operators (MPOs), and obtains the
occupied projector through density-matrix purification.

The package is intended for method development and controlled numerical
studies. It has small-system dense-reference tests, but it is not a
black-box production electronic-structure code: spectral bounds, MPO
truncation, convergence diagnostics, and sensitivity to initial seeds need
to be checked for every calculation.

## What is implemented

- Open one-dimensional chains with \(N=2^L\) sites.
- Open square lattices with \(L\) even and
  \(2^{L/2}\times2^{L/2}=2^L\) sites.
- Real scalar nearest-neighbour interaction \(U\), static potential `W`, and
  constant or position-dependent hopping `t`.
- Self-consistent Hartree and nearest-neighbour Fock fields in 1D and on the
  square lattice.
- Three zero-temperature purification modes:
  - `:sp2` — canonical trace-correcting second-order purification; the default.
  - `:palser_manolopoulos` — the historical canonical adaptive
    Palser–Manolopoulos/McWeeny path (`:adaptive_pm_mcweeny` remains an alias).
  - `:mcweeny_mu` — direct fixed-chemical-potential McWeeny purification.
- Observable and diagnostic reporting: density, bond order, particle number,
  energy components, Hermiticity, idempotency, and stationarity.
- Binary-carry square Hartree/Fock extraction, with TCI square Fock and
  experimental 1D alternatives retained for validation and comparison.

## Model and conventions

The one-body density matrix is

\[
\rho_{ij}=\langle c_j^\dagger c_i\rangle.
\]

Thus `rho[i,i]` is the occupation of site `i`. For an open nearest-neighbour
bond \(\langle i j\rangle\), the currently implemented real-exchange
Hartree–Fock functional is

\[
E = \operatorname{Tr}(H_0\rho) +
U\sum_{\langle ij\rangle}
\left[n_i n_j - \bigl(\operatorname{Re}\rho_{ij}\bigr)^2\right].
\]

Every physical bond is counted once. `H0` includes hopping and the external
potential `W`; there is no periodic wrapping or on-site interaction term.
The Fock field and energy intentionally use `real(rho[i,j])`, so the current
implementation is not a complex-exchange functional.

### Basis and coordinates

- Julia matrix index `i` corresponds to zero-based binary state `i - 1`.
  `MatrixChecker(M, sites, i, j, bra, ket)` evaluates `M[i,j]` in this basis.
- MPO sites are ordered from the most-significant binary bit to the least.
- In a 1D functional model, `t(x)` labels the open bond `x → x + 1` and
  `W(x)` labels site `x`, with zero-based `x`.
- In a square functional model, `t=(t_x,t_y)` and `W(x,y)` use zero-based
  coordinates. Coordinate bits are interleaved: `y` occupies the even bit
  positions and `x` the odd positions. `square_lattice_index`,
  `square_lattice_decoder`, and `square_neighbours` expose the package
  convention.

## Install

The project is tested with Julia 1.12.6.

```bash
git clone <repository-url>
cd MPO_Purification_HF
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Load it with:

```julia
using MPO_MeanField
```

The maintained solver workflow is CPU-first. CUDA remains in the project
environment for planned GPU work, but the current core package does not load a
graphics backend. A networked machine is normally required for the first
`Pkg.instantiate()`.

### Offline HPC environment

On an offline compute cluster, transfer a prepared Julia depot and point
`JULIA_DEPOT_PATH` to it. Do not use `Pkg.test()` on an offline node, because
it may attempt to initialise a registry. Instead run the maintained test
entry point directly:

```bash
module load julia/1.12.6
export JULIA_DEPOT_PATH=/gpfs/projects/<account>/cluster_depot
export JULIA_PKG_PRECOMPILE_AUTO=0
julia --startup-file=no --project=. test/runtests.jl
```

For the supplied Slurm job, see
[`benchmark/first_test/README.md`](benchmark/first_test/README.md).

## Quick start: non-interacting 1D calculation

All parameter fields are explicit. This is deliberate: numerical tolerances
and truncation limits are part of the calculation definition.

```julia
using MPO_MeanField

params = Parameters1D(
    L=2,                         # N = 2^2 = 4 sites
    t=-0.7,
    U=0.0,
    W=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
    S=nothing,                   # optional initial mean-field seed
    tci_tol=1e-10,
    itensors_tol=1e-12,
    itensors_maxdim=64,
    density=0.5,
    purification_steps=35,
    scf_mixing=0.5,
    scf_tol=0.1,
    scf_max_iterations=10,
)

sys = System(params)

# Arguments are (H_min, H_max), and must enclose every effective Hamiltonian
# encountered during the SCF calculation.
ok = run_scf!(sys, -5.0, 5.0;
    purification_method=:sp2,
    verify_spectral_bounds=true,  # exact CPU validation; only for N ≤ 16
    verbose=:all,
    record_energy=true,
)

ok || error("SCF did not converge: $(scf_diagnostics(sys).termination_reason)")
obs = observables(sys)
@show obs.particle_number obs.energy obs.idempotency_residual
```

`run_scf!` returns `true` only when its SCF residual history is stable.
Inspect `scf_diagnostics(sys)` even after a successful run.

## Quick start: open square lattice

For `ParametersSquare`, `L` must be even. `L=4` describes a \(4\times4\)
open lattice (16 sites), not a four-site lattice.

```julia
params = ParametersSquare(
    L=4,
    t=(-0.6, -0.35),             # (t_x, t_y)
    U=0.15,
    W=(x, y) -> 0.11x + 0.07y + 0.013x * y,
    S=nothing,
    tci_tol=1e-10,
    itensors_tol=1e-12,
    itensors_maxdim=64,
    density=0.5,
    purification_steps=35,
    scf_mixing=0.4,
    scf_tol=0.1,
    scf_max_iterations=40,
)

sys = System(params)
# Bounds include hopping, the known range of W, and a conservative bound on
# the square Hartree and real-exchange fields.
bounds = square_scf_spectral_bounds(
    params;
    potential_bounds=(0.0, 0.657), # extrema of W on this 4 x 4 example
    margin=0.5,
)
ok = run_scf!(sys, bounds...;
    purification_method=:sp2,
    verify_spectral_bounds=true,
    verbose=:all,
    record_energy=true,
)

diagnostics = scf_diagnostics(sys)
obs = observables_square(sys)
```

The corresponding tested weak-coupling reference, including diagnostics and
a density diagonal, is documented in
[`../docs/reference_cases/square_l4_weak_hf.md`](../docs/reference_cases/square_l4_weak_hf.md).

## Purification modes

All purification paths require a valid enclosing interval `(H_min, H_max)`.
For a small system, use `verify_spectral_bounds=true` to check it against the
dense spectrum at every SCF update. That validation is intentionally limited
to \(N\le16\) and must not be used as a production spectral solver.

For square SCF calculations with constant hopping, use
`square_scf_spectral_bounds` to construct a conservative interval from known
global coefficient bounds. For a functional `W`, supply its analytic range as
`potential_bounds=(W_min, W_max)`; for functional hopping, also supply
`hopping_abs_bounds=(max_abs_tx, max_abs_ty)`. The helper assumes a physical
density contraction and adds a worst-case `8abs(U)` interaction radius, so it
is deliberately conservative. Retain a positive `margin` for MPO truncation
uncertainty and document all supplied bounds in production metadata.

| Method | Ensemble / input | Use case |
| --- | --- | --- |
| `:sp2` | Canonical; target `round(N*density)` | Default method when the particle number is known. |
| `:palser_manolopoulos` | Canonical; target `round(N*density)` | Historical adaptive method; useful for comparisons. |
| `:mcweeny_mu` | Fixed chemical potential `chemical_potential=mu` | Grand-canonical zero-temperature projector \(\Theta(\mu-H)\). The trace is an output, not an imposed target. |

For `:mcweeny_mu`, `mu` must lie strictly inside the supplied spectral bounds.
If it lies in a gap, `Tr(rho)` should be close to an integer; at an eigenvalue
or a degeneracy, the zero-temperature particle number is ambiguous. The
direct McWeeny path now requires relative idempotency below `1e-6` before it
reports convergence.

For every method, check at least:

```julia
result = perform_purification(rho0, params;
    method=:sp2,
    spectral_bounds=(H_min, H_max),
    verbose=1,
)

@show result.converged result.termination_reason
@show result.trace result.trace_error
@show result.idempotency_residual result.hermiticity_residual
@show result.final_bond_dimension result.work
```

Do not continue an SCF calculation with an unconverged purification unless
you intentionally set `allow_unconverged_purification=true` and record why.

## Convergence and diagnostics

Purification diagnostics include the trace, idempotency residual,
Hermiticity residual, termination reason, iteration count, and bond-dimension
history. A canonical result needs both a trustworthy spectral interval and a
trace compatible with the intended filling; idempotency alone is not enough.

SCF convergence requires the Hartree-field, Fock-field, density, and
commutator residuals to remain below the configured threshold for consecutive
iterations (`stable_iterations=2` by default). The solver can stop with:

- `:converged` — stable SCF residual history;
- `:two_cycle_detected` — a density two-cycle rather than a fixed point;
- `:purification_failed` — purification stopped without converging;
- `:max_iterations` — SCF iteration budget exhausted.

`verbose=:all` prints updating terminal progress on a TTY; redirected logs
remain line-delimited. Use `record_energy=true` to retain the reported total
energy in the SCF history.

## Validation and tests

Run the test suite from the package directory:

```bash
julia --startup-file=no --project=. test/runtests.jl
```

The suite covers analytical and dense-reference limits, basis conventions,
Hamiltonian construction, field extraction, purification invariants, 1D and
square SCF references, energy double counting, and convergence diagnostics.
It also contains intentional difficult cases: exact Fermi-level degeneracy
must not be mistaken for a unique zero-temperature canonical projector.

The current suite records one known finite-gap square MPO-SP2 regression as
`@test_broken`. It is deliberately visible rather than treated as a validated
production result. Consult the test output and
[`../docs/plans/remaining_implementation_priorities.md`](../docs/plans/remaining_implementation_priorities.md)
before selecting SP2 settings for a new large square calculation.

## Benchmarks

`benchmark/extraction_scaling.jl` compares 1D Hartree/Fock extraction paths
over `L=4:2:14` by default and writes machine-readable samples, summaries,
bond dimensions, and process-memory information. It is intended for a
separate CPU cluster run, not routine development.

```bash
mkdir -p benchmark_results/manual
bash benchmark/run_extraction_scaling.sh benchmark_results/manual
```

The Slurm template is
[`benchmark/extraction_scaling.slurm`](benchmark/extraction_scaling.slurm).
Set its account, QoS, module name, and `JULIA_DEPOT_PATH` appropriately before
submitting it. Read the generated `README.md`, `samples.csv`, `summary.csv`,
and `process_time.txt` together: peak RSS is process-wide, while the CSVs hold
per-method timing and allocation data.

## Current limitations

- Open boundaries only; no periodic-boundary implementation.
- Spinless, real scalar nearest-neighbour `U` only; no multiorbital, spin,
  long-range, or spatially varying interaction model.
- The present exchange functional keeps only the real nearest-neighbour bond
  order; complex exchange physics is not represented.
- `:sp2` can stagnate under MPO truncation in some finite-gap square cases.
  Check `PurificationResult` and do not rely on a convergence flag alone.
- Exact spectral-bound validation is a small-system diagnostic, not a scalable
  bound estimator. Production calculations need conservative externally
  justified bounds.
- GPU movement hooks exist in `run_scf!`, but the maintained validation and
  benchmark workflow is CPU-first. Treat GPU use as experimental until it has
  its own validation run.

## Repository map

```text
src/
  core/            Parameter types, system state, and operators
  hamiltonians/    Open-chain and open-square MPO Hamiltonian construction
  purification/    Initial scaling, SP2, direct McWeeny, and PM backends
  tci/             TCI and tensorial mean-field extraction
  hf/              Self-consistent-field driver and diagnostics
  utils/           Quantics helpers, progress, energies, and observables
test/              Small-system, dense-reference, and regression tests
benchmark/         CPU scaling and cluster-job scripts
```

For the active scientific and implementation roadmap, see
[`../docs/plans/remaining_implementation_priorities.md`](../docs/plans/remaining_implementation_priorities.md).
