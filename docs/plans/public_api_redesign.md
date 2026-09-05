# Public API redesign plan: model, representation, solver, result

## Scope and non-goals

This is an architecture and migration plan for the public interface of
`MPO_MeanField`. It does **not** change Hartree--Fock conventions, QTT
indexing, MPO truncation, SP2, KPM, or any existing scientific kernel.

The goal is to make a calculation readable in physical terms, while retaining
every numerical control required for a reproducible tensor-network study:

\[
\text{Model} \longrightarrow \text{Representation} \longrightarrow
\text{Solver} \longrightarrow \text{Result}.

## Migration status

The compatibility migration is complete enough for the first public production
MPO campaign: an equal open square, CPU/CUDA runtime choice, explicit
spectral enclosure, reproducible QTT settings, SP2 and Fock controls, and
normalized result output are all expressed through public configuration.

In particular, `SCFSettings` now owns the operational SP2 idempotency and
relative trace tolerances, square Fock backend, energy-history policy,
stability/two-cycle policy, and safeguarded linear/Pulay mixing controls. The
first public (64\times64) separable
quasiperiodic CUDA campaign is
`runs/2d/separable_aubry_andre_lside6/campaign_public.jl` with its
site-specific launcher `submit_public_mpo_cuda.slurm`.

This is a compatibility layer over the proven legacy equal-square MPO kernel;
it is not a claim that rectangular MPO/dense support, alternative QTT
encodings, automatic spectral enclosures for arbitrary functions, or all
legacy diagnostics have been migrated. Those remain explicit follow-up work.
\]

The current low-level API is scientifically appropriate for expert use. The
redesign is not intended to conceal numerical settings or turn the package
into an unvalidated black box. It instead removes accidental coupling between
physical parameters, QTT representation, solver controls, runtime selection,
and output layout.

## 1. Current API map

### Current public package API

The package exports the following user-facing layers from
`src/MPO_MeanField.jl`:

| Layer | Existing objects/functions | Current responsibility |
|---|---|---|
| Model plus numerics | `Parameters1D`, `ParametersSquare` | Mixes geometry, hopping, interaction, potential, seed, filling, QTT/TCI tolerance, MPO cutoff/cap, and SCF controls. |
| State construction | `System(params)` | Validates parameters; makes qubit indices, translations, `H0`, initial `VH`, `VF`, and `rho`; precomputes basis states. |
| Solver | `run_scf!(sys, H_min, H_max; ...)` | Runs purification, mean-field extraction, linear mixing, diagnostics, and stopping logic. |
| Purification | `construct_rho_0`, `perform_purification`, SP2/PM/McWeeny backends | Scales the Hamiltonian and produces a density-projector MPO. |
| Representation utilities | `square_lattice_index`, `square_lattice_decoder`, `MatrixChecker`, QTT/TCI construction | Encodes a geometry in a binary MPO basis and accesses entries. |
| Observables | `observables`, `observables_1d`, `observables_square` | Returns densities, bonds, energies, particle number, and residuals from a `System`. |

The core numerical flow is:

```text
Parameters* → System → H_eff = H0 + VH + VF
                         ↓
                  construct_rho_0
                         ↓
                  perform_purification
                         ↓
                 extract Hartree/Fock fields
                         ↓
                   SCF residuals + mixing
                         ↓
                 observables + SCFDiagnostics
```

The implementation is principally in `src/core/system.jl`,
`src/hamiltonians/mpo_construction.jl`, `src/purification/`,
`src/hf/self_consistent.jl`, and `src/utils/observables.jl`.

### Existing geometry and QTT rules

`Parameters1D(L=...)` means a binary chain of length

\[
N = 2^L.
\]

`ParametersSquare(L=...)` requires even `L` and means an open square with

\[
N_x=N_y=2^{L/2},\qquad N=2^L.
\]

The square convention is fixed interleaved coordinate bits: integer bit
positions `0,2,4,...` encode `y`, and `1,3,5,...` encode `x`. The MPO index
order itself is most-significant binary bit first. These conventions are used
by Hamiltonian construction, binary-carry fields, local measurement, and
dense-reference scripts; they cannot change during this API migration.

### Existing runners and output contracts

There is a sound common contract for explicit 1D/square campaign runners:

```text
input.jl, selection.toml, metadata.toml,
scf_history.csv, observables.toml, site_density.csv, bond_order.csv
```

However, the repository also contains numerous diagnostic and KPM runners
with different positional CLI signatures and output names such as
`summary.toml`, `iterations.csv`, `sp2_progress.txt`, and task-specific CSV
files. This is the primary usability problem: the same physical calculation
is described differently in an interactive session, a regular campaign, a
fixed-H validation, and a KPM workflow.

## 2. Parameter classification and propagation

The table gives the target ownership of every currently user-selectable
parameter. “Current propagation” describes existing behavior to preserve in
the compatibility layer.

| Current parameter/control | Target layer | Current propagation | Notes |
|---|---|---|---|
| `Parameters1D.L` | Model | Sites, translations, `N`, all local loops | Replace publicly with `size=N`; derive QTT bits. |
| `ParametersSquare.L` | Model plus derived representation metadata | Equal square side, interleaved indexing, `N`, all square loops | Replace publicly with `size=(Nx,Ny)`; the model accepts rectangles, while legacy kernels remain equal-square until explicitly generalized. |
| `t` | Model | `build_H0`; dense and KPM runners construct equivalent hopping separately | Physical bond hopping. 1D accepts scalar/function; square accepts `(tx,ty)`, each scalar/function. |
| `U` | Model | Hartree/Fock extraction and energy | Current physics is real scalar nearest-neighbour repulsion/attraction. |
| `W` | Model | TCI diagonal MPO / dense diagonal | Static external potential. |
| `S` | Model | Initial Hartree field / dense seed | Initial symmetry-breaking field; must remain distinct from `W`. |
| `density` | Model | `Ne=round(N*density)` in canonical methods | A physical filling, although rounding convention is numerical provenance. |
| boundary condition | Model | Implicitly open in all current core code | Make explicit and initially validate only `:open`. |
| `tci_tol` | Representation | Function-to-QTT diagonal MPO construction | A representation/interpolation error control, not model physics. |
| `itensors_tol` | Representation | Every MPO `apply`, sum, and compression | Rename public field to `cutoff`; record full effective value. |
| `itensors_maxdim` | Representation | Every MPO `apply`, sum, and compression | Rename public field to `maxdim`. |
| QTT bit ordering | Representation | Implicit fixed 1D binary / square interleaved convention | Expose as declared metadata; do not expose unsupported alternate encodings. |
| `purification_steps` | Solver | Maximum inner purification iterations | Keep explicit. |
| `purification_method` | Solver | `run_scf!` → `perform_purification` | SP2, Palser--Manolopoulos, and fixed-\(\mu\) McWeeny are core methods. |
| `scf_mixing` | Solver | Linear mixing of Hartree and Fock MPOs | Current core MPO SCF uses linear mixing only. |
| `scf_tol` | Solver | Converted internally to `scf_tol/100` residual threshold | Preserve this exact behavior in the first wrapper; document the resulting residual tolerance. |
| `scf_max_iterations` | Solver | Outer SCF loop bound | Keep explicit. |
| `stable_iterations`, two-cycle controls | Solver | `run_scf!` stopping policy | Advanced settings, but meaningful and reproducible. |
| SP2 idempotency/trace tolerance | Solver | `perform_purification_sp2` stopping policy | Solver accuracy settings, separate from MPO cutoff. |
| `chemical_potential`, McWeeny trace options/form | Solver | Fixed-\(\mu\) purification only | An ensemble/algorithm setting; model still owns filling for canonical runs. |
| supplied `(H_min,H_max)` | Solver bound policy | Initial scaling at each SCF iteration | Replace raw positional arguments by a policy object. |
| exact bound validation, safety margin | Solver validation / bound policy | Optional small-system check in `construct_rho_0` | Must remain visible in results. |
| `square_fock_method` | Solver implementation setting | Chooses `:binary_carry` or `:tci` extraction | Same intended functional; numerical backend must be recorded. |
| garbage-collection controls | Runtime | Purification cleanup cadence | Not physical and normally advanced. |
| CPU/CUDA backend, precision, device setup | Runtime | Runner-specific callbacks and separate KPM implementations | Must become an early, explicit runtime choice. |
| verbosity, progress IO, phase callback | Runtime/output | Progress, profiling, diagnostics | Separate from numerical settings. |
| result root, campaign name, Slurm info, output detail | Output/provenance | Campaign runners only | Never change model or solver behavior. |

## 3. Proposed public types

The following types are public because they correspond to choices a user must
understand, compare, serialize, or cite. Their internals should be immutable
and printable. Construction should validate immediately and return helpful
errors.

```julia
abstract type AbstractPhysicalModel end

ChainModel(; size, hopping, interaction, potential=nothing, seed=nothing,
             filling, boundary=:open)

SquareModel(; size::Tuple{Int,Int}, hopping, interaction,
              potential=nothing, seed=nothing, filling,
              boundary=:open)

QTTSettings(; encoding=:interleaved, tci_tol=1e-10,
              cutoff=1e-10, maxdim=256)

SCFSettings(; purification=:sp2, mixing=0.5, tolerance=0.1,
              maxiter=30, purification_maxiter=50,
              bounds=AutoBounds(), stable_iterations=2,
              detect_two_cycles=true, square_fock_method=:binary_carry,
              sp2_idempotency_tolerance=1e-3)

RuntimeSettings(; backend=:cpu, device_scalar_type=Float64,
                 threads=1, progress=:batch)
```

### `ChainModel` and `SquareModel`

These belong in the public API because physical size, fields, interaction,
filling, and boundary conditions define the scientific system independently
of MPO, dense, or KPM treatment.

The public geometry accepts all open, power-of-two rectangles. Solver support
is deliberately narrower and must be checked independently:

| Public model | Initial support | Rejection rule |
|---|---|---|
| `ChainModel(size=N)` | Open chain with `N=2^L` | Reject non-power-of-two `N`. |
| `SquareModel(size=(Nx,Ny))` | Open square or rectangle with `Nx=2^ell_x`, `Ny=2^ell_y` | Reject non-powers of two in the constructor; legacy MPO/dense conversion rejects `Nx != Ny` until those kernels are generalized. |

This avoids conflating physical geometry with current solver availability.
`SquareModel(size=(128,64))` is a valid model specification and exposes its
derived QTT bits. A backend preflight must nevertheless reject a requested
MPO/dense calculation until the Hamiltonian, coordinates, local observables,
dense reference, and MPO field extractors support unequal bit counts end to
end. The existing rectangular KPM runner is useful implementation evidence,
not yet a generic public backend.

### `QTTSettings`

This belongs in the public API because cutoff, cap, and interpolation error
are a controlled approximation whose values affect the result. It should not
pretend there are interchangeable encodings:

- Chain models initially use the current binary ordering; `encoding=:binary`
  is metadata, not a new behavior.
- Square models initially accept only `encoding=:interleaved`.
- Unsupported choices such as a blocked `x`-then-`y` encoding must fail in
  preflight instead of falling back silently.

Existing `tci_tol`, `itensors_tol`, and `itensors_maxdim` map directly to
`QTTSettings.tci_tol`, `.cutoff`, and `.maxdim` in the first compatibility
layer.

### `SCFSettings` and bound policies

This belongs in the public API because purification, stopping criteria,
mixing, and spectral enclosure define the numerical algorithm. It must not
live in the model.

Use a typed bound policy rather than raw positional `H_min,H_max`:

```julia
AutoBounds(; margin=0.5,
             potential_bounds=nothing,
             hopping_abs_bounds=nothing,
             validate_small_system=false)

ExplicitBounds(lower, upper; margin=0.0,
               validate_small_system=false)
```

`AutoBounds` invokes the existing `square_scf_spectral_bounds` only when its
assumptions can be satisfied. It may automatically cover constant hopping,
no potential, and known built-in fields. It must not attempt to infer extrema
of an arbitrary Julia function by sampling it. For a custom function,
preflight requests analytic envelopes via `potential_bounds` and, for
functional hopping, `hopping_abs_bounds`; otherwise it refuses the run.

The fully expanded numerical interval and source (`automatic`, `explicit`,
or `exact_small_system_validated`) must be written to provenance. The first
implementation continues to pass a fixed enclosing pair to existing
`run_scf!`; dynamic spectral-bound updates are explicitly out of scope.

### `RuntimeSettings`

This belongs in the public API because backend selection determines
availability, data movement, precision, memory limits, and expected
validation. It should be requested at the beginning of a command or
`solve` call, not hidden in Slurm glue.

The first wrapper may expose `backend=:cpu` and `backend=:cuda`, but must
report which path is actually used. Current MPO SCF uses GPU movement hooks
but returns purified density to host for mean-field extraction; it is not an
end-to-end GPU mean-field implementation. The KPM GPU solvers are separate,
more device-oriented implementations. Preflight must state this distinction.

## 4. Minimal and expert APIs

### Minimal valid API

```julia
using MPO_MeanField

model = SquareModel(
    size=(32, 32),
    hopping=aubry_andre_hopping(V2=0.5),
    interaction=0.5,
    filling=0.5,
)

result = solve(model)
```

`solve(model)` uses documented, inspectable defaults, equivalent to calling
`solve(model; representation=default_representation(model),
solver=default_solver(model), runtime=RuntimeSettings())`. It must print or
return the expanded configurations; it must never hide them. For a generic
functional hopping, the minimal call may legitimately fail at preflight until
the user supplies a safe bound envelope. This is preferable to manufacturing
an unsafe spectral interval.

The helper `aubry_andre_hopping` must be a public built-in model field only
after its phase and bond convention are documented and tested. Before that,
the example should use a constant hopping or an explicit function plus an
`AutoBounds` envelope.

### Expert API

```julia
model = SquareModel(
    size=(64, 64),
    hopping=(tx, ty), interaction=1.0, potential=nothing,
    seed=checkerboard_seed(0.1), filling=0.5, boundary=:open,
)

representation = QTTSettings(
    encoding=:interleaved, tci_tol=1e-10, cutoff=1e-10, maxdim=512,
)

solver = SCFSettings(
    purification=:sp2, mixing=0.5, tolerance=0.1, maxiter=30,
    purification_maxiter=50,
    bounds=AutoBounds(
        hopping_abs_bounds=(1.5, 1.5), margin=0.5,
        validate_small_system=false,
    ),
    sp2_idempotency_tolerance=2e-4,
)

runtime = RuntimeSettings(backend=:cuda, device_scalar_type=Float64)
preflight_report = preflight(model; representation, solver, runtime)
result = solve(model; representation, solver, runtime)
```

Every proposed keyword is public only because it affects the physical system,
controlled approximation, algorithm, or reproducibility. Low-level ITensor
site construction, precomputed basis states, GC callbacks, `to_gpu`/`to_cpu`
functions, and raw MPO objects remain internal implementation details.

## 5. Physical-size to QTT-size rules

| Model | User supplies | Derived internal bits | Initial support |
|---|---|---|---|
| Chain | `size=N=2^ell` | `L_qtt=ell` | Yes |
| Square | `size=(2^ell,2^ell)` | `Lx=Ly=ell`, `L_qtt=2ell` | Public model; legacy MPO/dense conversion available |
| Rectangle | `size=(2^ell_x,2^ell_y)` | `L_qtt=ell_x+ell_y`; low bits are interleaved, then unmatched coordinate bits are appended in documented order | Public model; no generic MPO/dense backend yet |

For the existing square implementation:

\[
(N_x,N_y)=(64,64) \Rightarrow (L_x,L_y)=(6,6),\quad L_{\rm QTT}=12.
\]

The new model constructors should calculate these values and store them in a
derived `GeometryInfo`, returned by `preflight` and recorded in the result.
The constructors should reject `size=(48,48)` rather than round it, pad it,
or choose a representation without user approval.

## 6. Preflight design

```julia
report = preflight(model; representation, solver, runtime)
show(report)
```

`preflight` is public because it makes expensive numerical assumptions
inspectable before execution. It performs no SCF iteration and does not
alter result directories. At minimum it reports:

```text
Geometry                 square, open
Physical lattice         64 × 64
Physical sites           4096

QTT bits / x, y          6, 6
Total QTT levels         12
Encoding                 interleaved (current supported square encoding)

Filling                  0.5
Target particles         2048 (round(N × filling))

Purification             canonical SP2
MPO cutoff               1e-8
Maximum bond dimension   512

Spectral bounds          [-8.5, 12.5]
Bound source             auto: square envelope + margin
Backend                  CUDA / Float64
```

Warnings must distinguish hard errors from non-fatal risks:

| Condition | Preflight action |
|---|---|
| Non-power-of-two physical size | Error for initial QTT API. |
| Rectangular geometry requested with legacy MPO/dense solver | Error before construction/allocation; explain that the model is valid but this backend is not yet generalized. |
| Unknown envelope for arbitrary `W` or hopping function | Error unless explicit bounds are supplied. |
| `backend=:cuda` but CUDA unavailable | Error before model construction/allocation. |
| GPU MPO-SCF request | Warning that mean-field extraction currently returns to CPU. |
| Square SP2 with settings outside existing validation coverage | Warning, including the selected cap/cutoff and recommendation for fixed-H/dense calibration. |
| `verify_small_system=true` with large N | Error or explicit disable with explanation. |
| Estimated dense matrix too large | Warning/error based on dense backend policy. |
| `maxdim` or requested method likely expensive | Warning based on documented heuristic only; never a claimed memory guarantee. |

## 7. Solver-independent result design

The common result must represent the *physical answer* and numerical status,
not require every backend to materialize identical internal objects.

```julia
struct SolveResult{M,R,S}
    model::M                         # expanded, immutable physical model
    representation::R                # expanded representation settings
    solver::S                        # expanded solver settings and bounds used
    runtime::RuntimeProvenance
    status::SolveStatus              # converged, termination reason, warnings
    observables::PhysicalObservables
    history::Vector{SCFHistoryRecord}
    diagnostics::AbstractSolverDiagnostics
    backend_state                    # optional expert-only state; may be nothing
end
```

`PhysicalObservables` should provide consistent names where available:

```julia
density
bond_order
particle_number
energy = (kinetic, hartree, fock, interaction, total)
hermiticity_residual
idempotency_residual
stationarity_residual
```

Values unavailable to a solver are `missing`, not zero and not silently
estimated. Large solvers may offer local observables through arrays, sampled
sets, or storage-backed handles; provenance must state coverage. A KPM
result must not advertise a full density if only selected entries were
computed.

`SCFHistoryRecord` has common fields with `missing` for inapplicable values:
iteration, trace, trace error, idempotency defect, field/density/bond
residuals, stationarity residual, total energy, and effective numerical
controls. SP2, KPM, and dense-specific records remain in
`diagnostics`:

| Backend | Solver-specific diagnostics |
|---|---|
| Dense | eigenspectrum/eigensolver timing, exact stationarity |
| MPO-SP2 | branch history, purification trace/idempotency history, bond-dimension and truncation history |
| KPM | moments, kernel, probes, blocking, chemical-potential solve, stochastic/nested-probe audit |

This type belongs in the public API because it is the comparison boundary for
physical results and the foundation for a stable output schema. Raw `System`,
MPOs, CUDA arrays, and Chebyshev work buffers are not portable public result
fields.

## 8. Dense, MPO, and KPM compatibility

| Capability | Dense campaign | MPO-SP2 core | Square KPM | Rectangular KPM scripts |
|---|---:|---:|---:|---:|
| 1D open chain | Yes | Yes | No generic implementation | No |
| Equal open square | Yes | Yes | Yes | N/A |
| Arbitrary rectangular lattice | No | No | No generic implementation | A specific hard-coded/runner-level path exists |
| Scalar/function hopping | Yes for current parameter types | Yes | Uses `ParametersSquare` in square runner | Separate hard-coded model data |
| Static potential/seed | Yes | Yes | Yes in square runner | Runner-specific |
| Real scalar nearest-neighbour U | Yes | Yes | Yes | Runner-specific |
| Canonical filling | Exact occupied states | SP2/PM canonical | KPM chemical-potential solve | KPM chemical-potential solve |
| Fixed chemical potential | Not a campaign public option | McWeeny-\(\mu\) | KPM implementation-specific | KPM implementation-specific |
| Full density/bonds | Yes, practical only at small N | Yes, measurement can be costly | Yes in current square KPM SCF | Depends on runner mode |
| GPU | Not the maintained dense path | Experimental movement hooks; extraction remains host-side | CUDA supported and validated by dedicated workflows | CUDA supported by dedicated workflows |

Therefore a single call form is realistic, but not a false common physics
abstraction:

```julia
solve(model; method=:dense, ...)
solve(model; method=:mpo, ...)
solve(model; method=:kpm, ...)
```

Initial implementation should dispatch only combinations in the table that
are already supported. A `MethodCompatibilityError` must explain unsupported
combinations before a calculation begins. Dense and MPO wrappers are
straightforward because they already consume the same `Parameters1D` or
`ParametersSquare`. KPM requires a substantial adapter because its runners
currently own independent configuration, local-observable kernels, chemical
potential logic, GPU setup, and output writing.

## 9. Campaign integration

Use the same objects interactively and on Slurm. A campaign should be a Julia
source file, not pure TOML, because hopping, potential, and seed are often
functions.

```julia
campaign = CampaignSpec(
    name="separable_aa_64x64_validation",
    cases=[
        CaseSpec(
            label="V2_0p5_U_1",
            model=SquareModel(...),
            representation=QTTSettings(...),
            solver=SCFSettings(...),
            validation=ValidationSettings(dense_reference=true,
                                          fixed_h_sp2=true),
        ),
    ],
    runtime=RuntimeSettings(backend=:cuda),
    output=OutputSettings(root=ENV["MPO_RESULTS_ROOT"]),
)
```

The generic runner is conceptually:

```bash
julia --project=. bin/mpohf.jl preflight campaign.jl --task 1 --method mpo --backend cuda
julia --project=. bin/mpohf.jl run       campaign.jl --task 1 --method mpo --backend cuda
```

Slurm passes only the task index, method, and explicit backend. It does not
encode physics or numerics as a fragile list of positional Julia arguments.
Expert fixed-H and compression studies become named campaign stages, not
unrelated executable names.

All stages write schema-versioned output:

```text
run.toml                 # expanded model/representation/solver/runtime
status.toml              # convergence status and warnings
observables.toml         # common physical scalars
convergence.csv          # normalized outer/inner trajectory
density.csv              # coverage recorded in run.toml
bond_order.csv           # coverage and convention recorded in run.toml
timings.toml
diagnostics/             # method/stage-specific detail
```

Existing immutable `input.jl`, Slurm metadata, and refusal to overwrite
results are retained. `run.toml` complements rather than replaces `input.jl`:
the former records the expanded canonical configuration; the latter preserves
the executable function definitions.

## 10. Backward compatibility

The first implementation wraps the existing validated path:

```text
SquareModel + QTTSettings + SCFSettings
                 ↓ conversion
ParametersSquare
                 ↓
System(params)
                 ↓
run_scf!(sys, lower, upper; existing keywords)
                 ↓
observables_square(sys), scf_diagnostics(sys)
                 ↓ conversion
SolveResult
```

Similarly, `ChainModel` maps to `Parameters1D`. Existing public exports,
campaign files, and runners remain supported unchanged for at least one
release cycle. They should be documented as “legacy expert API”, not removed.

The following remain internal in the new design:

- `System` layout and mutable MPO state;
- `build_H0`, translation MPOs, precomputed qubit states;
- `MatrixChecker` and QTT coordinate implementation;
- `construct_rho_0` and individual purification work arrays;
- `to_gpu`, `to_cpu`, phase callbacks, GC callbacks;
- raw KPM blocks/probe arrays.

They may stay exported temporarily for backwards compatibility, but new
documentation should not lead users there for standard calculations.

## 11. Rectangular MPO enablement

Rectangular geometry is now a valid public model specification, but it is not
yet an MPO/dense solver capability. Enabling it must be a separate
scientific-kernel project, not an incidental consequence of the API redesign.
It must preserve the existing equal-square indexing exactly.

### 11.1 Canonical rectangular bit convention

For

\[
N_x=2^{L_x},\qquad N_y=2^{L_y},\qquad L_{\rm QTT}=L_x+L_y,
\]

use the convention already exercised by the rectangular KPM workflow:

1. interleave the common low coordinate bits in integer-bit order
   \((y_0,x_0,y_1,x_1,\ldots)\);
2. append any remaining higher bits of the longer coordinate in increasing
   significance order;
3. retain the existing MPO tensor order, from most-significant QTT bit to
   least-significant QTT bit.

When \(L_x=L_y\), this reduces identically to the existing square convention.
The convention must be documented in code and result provenance; no implicit
padding, bit reordering, or alternate encoding is permitted.

### 11.2 Implementation order and acceptance tests

| Rectangular stage | Change | Required verification before proceeding |
|---|---|---|
| R1: geometry helpers | Add rectangular `index`, `decoder`, neighbour and bond helpers parameterized by `(Lx,Ly)`; retain existing square helpers as compatibility wrappers. | Exhaustively prove `decode(index(x,y))=(x,y)` and correct open-boundary neighbours/bond count on `2×2`, `4×2`, `2×4`, `8×4`, and `4×8`. Equality with current square helpers is exact for `4×4` and `8×8`. |
| R2: hopping MPO | Generalize `build_translation_square` (or replace internally by a geometry-parameterized builder) to independent x/y active-bit sets and coordinate-specific carry boundaries. | Materialize the MPO densely for `4×2`, `2×4`, `8×4`, and `4×8`; compare all four translations and both constant/function-hopping Hamiltonians elementwise with independently written dense references. Verify Hermiticity and open boundaries. |
| R3: static diagonal terms and state construction | Generalize `build_W`, seed construction, `System`, and parameter validation without changing equal-square behavior. | Square regression tests remain identical at their current tolerances. Rectangular dense Hamiltonian includes every onsite and nearest-neighbour term exactly once. |
| R4: Hartree and Fock extraction | Generalize binary-carry shifts, adjacency Hartree MPO, Fock-band extraction, and local observables using the R1 coordinate map. | On small rectangles, compare every Hartree site, horizontal/vertical Fock bond, mean-field energy, trace, and Hermiticity against direct dense measurement. Test both orientations with unequal bit counts. |
| R5: purification and SCF | Permit rectangular `Parameters`/compatibility conversion, spectral bounds, SCF diagnostics, outputs, and dense reference solver. | Non-interacting MPO--dense agreement and interacting fixed-point MPO--dense agreement on at least `4×2` and `4×4`; test both x- and y-oriented modulation. |
| R6: public availability | Allow `solve(...; method=:mpo/:dense)` and campaign preflight for rectangles; add a KPM adapter only where its model mapping is independently verified. | Public rectangular example emits the common output schema and records the coordinate encoding, backend, tolerances, bounds, and validation status. |

R4 is the first stage that changes the mean-field implementation. It must not
begin until R1--R3 establish an exact non-interacting rectangular Hamiltonian.
The field-extraction tests are particularly important because the existing
binary-carry code assumes alternating x/y bits at every level; this is no
longer true after one coordinate runs out of bits.

## 12. API migration stages and required tests

No stage changes physical equations, QTT order, purification recurrence, or
truncation policy. Each should be a small reviewable commit.

| Stage | Change | Acceptance criterion |
|---|---|---|
| 1 | Add immutable model/settings types and conversion functions to legacy `Parameters*`; do not add `solve`. | Conversion reproduces every legacy parameter byte-for-byte; 1D and square constructor validation tests pass. |
| 2 | Add geometry/QTT conversion and `preflight`. | Power-of-two rejection, rectangular geometry/QTT metadata, backend capability errors, particle count, and supplied/automatic bound reports have deterministic tests. |
| 3 | Add `SolveResult` and a `solve(...; method=:mpo)` compatibility wrapper. | Existing small 1D and square reference tests reproduce density, energy, trace, residuals, termination reason, and SCF history. |
| 4 | Add common writer/reader schema around the MPO wrapper; retain old files. | `write_result`/`read_result` round-trip density, bonds, scalar observables, and SCF history for a small reference; public outputs use format-versioned `input.toml`, `metadata.toml`, `observables.toml`, `scf_history.csv`, `site_density.csv`, and `bond_order.csv`. Function-valued fields are descriptive in `input.toml`, not executable source. |
| 5 | Add `method=:dense` wrapper and dense compatibility checks. | Complete: the public dense wrapper directly diagonalizes the real open-boundary nearest-neighbour (H_{\rm eff}), occupies the canonical lowest (N_e) orbitals, and uses the legacy linear-mixing/stability policy. Deterministic 1D and interacting (4\times4) MPO--dense comparisons agree at the existing small-system tolerances. Dense results declare spectral bounds `not_used`; supplied bounds are rejected in preflight. |
| 6 | Add `CampaignSpec` and generic `preflight`/`run` CLI, then migrate one 1D and one equal-square campaign. | Complete: `bin/mpohf.jl` selects one `CaseSpec` by task and dispatches the validated public MPO/dense methods. Public 1D and equal-square separable-Aubry--André campaign sources live alongside, rather than replacing, their legacy expert campaigns. Runs are non-overwriting, write the common public schema, and copy immutable `campaign.jl` provenance. |
| 7 | Add explicit runtime/backend preflight and standardized timing/provenance. | Complete: `RuntimeSettings` and `runtime_preflight` report requested versus active backend, CUDA device information when usable, and the actual `cuda_purification_host_fields` transfer path. CUDA is lazy-loaded and unavailable paths fail before `System` construction. Public output records runtime provenance and whole-`solve` elapsed time. |
| 8 | Add KPM adapter only for already validated equal-square cases. | Complete: `method=:kpm` exposes Jackson moments, deterministic Hadamard probes, a canonical trace/chemical-potential solve, and a separate fixed-field audit through `KPMSettings`. It writes the common result schema and has a small equal-square KPM--dense regression. The public adapter currently uses linear SCF mixing; Pulay remains an expert-runner workflow pending its own validation/migration. |
| 9 | Generalize the public KPM graph adapter to open power-of-two rectangles. | Complete: `SquareModel(size=(Nx,Ny))` now maps directly to a row-major sparse nearest-neighbour graph, while a coordinate-interleaved Hadamard row code preserves the nested two-dimensional probe hierarchy. Small `4×2`, `2×4`, and `8×4` graph checks plus a `4×2` KPM--direct-diagonalization regression verify this adapter. This does not generalize legacy MPO or dense rectangle kernels. |
| 10 | Expose a geometry-agnostic sparse-graph KPM model. | Complete: `GraphModel` accepts a real symmetric sparse hopping matrix, explicit onsite/seed vectors, a global nearest-neighbour interaction, and optional Hadamard row codes. It is deliberately KPM-only; preflight rejects MPO/dense use. A four-site non-chain graph verifies sparse assembly, explicit probe-code persistence, KPM solve/output, and rejection of asymmetric graphs. |

Rectangular MPO/dense support follows R1--R6 above rather than being folded
into the API migration. Rectangular KPM adaptation, alternate encodings,
DIIS/Pulay for MPO, dynamic bounds, checkpointing, and any new physical model
are separate future projects.

## 13. Recommended implementation sequence

Implement stages 1--3 first. They produce immediate interactive clarity with
no campaign rewrite and establish the model/representation/solver/result
boundary. Pause there to use the new wrapper on existing small references.

Only after exact behavioral equivalence is demonstrated should stages 4--6
unify campaigns and outputs. CUDA standardization is stage 7 because backend
policy must be defined around actual data movement rather than a superficial
`backend=:cuda` keyword. KPM is last: it can share physical model concepts and
results, but it should not be forced into an MPO-shaped implementation before
its current separate kernels and diagnostics are deliberately adapted.

After API stage 3 has reproduced existing equal-square references, start
R1--R3 as one non-interacting rectangular-Hamiltonian milestone. Do not expose
rectangular MPO solving after that milestone: R4 and R5 are required to
validate the Hartree--Fock functional and self-consistency. This ordering
keeps the equal-square conference workflow stable while rectangular support is
demonstrated rather than assumed.

### Stages 6--7 public commands

The two migrated sources are
`runs/1d/aubry_andre_nn_l10/campaign_public.jl` and
`runs/2d/separable_aubry_andre_lside5/campaign_public.jl`.

```bash
julia --project=. bin/mpohf.jl preflight \
  runs/2d/separable_aubry_andre_lside5/campaign_public.jl \
  --task 2 --method mpo --backend cpu

julia --project=. bin/mpohf.jl run \
  runs/1d/aubry_andre_nn_l10/campaign_public.jl \
  --task 1 --method dense --backend cpu \
  --output-root "$MPO_RESULTS_ROOT"
```

For an allocated GPU node, request CUDA explicitly. The preflight reports the
hybrid execution path before it allocates an MPO; it does not claim that
Hartree/Fock extraction is on the device:

```bash
julia --project=. bin/mpohf.jl preflight \
  runs/2d/separable_aubry_andre_lside5/campaign_public.jl \
  --task 2 --method mpo --backend cuda
```

The corresponding generic Slurm launcher carries only campaign location,
method, and one-based array task; it does not replicate physics parameters:

```bash
export MPOHF_CAMPAIGN=runs/2d/separable_aubry_andre_lside5/campaign_public.jl
export MPOHF_METHOD=mpo
sbatch --array=1-3 --export=ALL,MPOHF_CAMPAIGN,MPOHF_METHOD,MPO_RESULTS_ROOT \
  runs/submit_public_campaign.slurm
```

`--task` is one-based, matching Slurm array task identifiers. A public
campaign case may retain MPO spectral bounds for `method=:mpo`; the same case
can be run with `method=:dense`, which ignores no hidden configuration and
therefore omits those bounds explicitly.
