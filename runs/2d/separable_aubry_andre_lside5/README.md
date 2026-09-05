# Separable 2D Aubry--Andre HF, side level 5

This campaign is the square-lattice counterpart of the validated 1D
quasiperiodic reference cases. `side_level=5` means a `32 x 32` open square:
`ParametersSquare(L=10)`, `N=1024`, and `N_e=512` at half filling.

The hopping is separable by direction:

\[
t_x(x,y)=-1-V_2\cos[2\pi\tau(x+1/2)],\qquad
t_y(x,y)=-1-V_2\cos[2\pi\tau(y+1/2)],
\]

where \(\tau=\sqrt2-5/6\). There is no external onsite potential. The
temporary seed `S(x,y)=0.1(-1)^{x+y}` appears only in the first SCF
Hamiltonian and selects a checkerboard branch; it is replaced by the
self-consistent Hartree field.

| Array task | `V2` | `U` | Conservative SP2 bounds |
| --- | ---: | ---: | --- |
| 1 | 0.0 | 1.0 | `(-12.5, 12.5)` |
| 2 | 0.5 | 1.0 | `(-14.5, 14.5)` |
| 3 | 2.0 | 0.2 | `(-14.1, 14.1)` |

Submit the MPO calculation from the package root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch runs/2d/separable_aubry_andre_lside5/submit.slurm
```

Then submit the independent dense one-particle HF reference:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_dense_hf.slurm
```

The dense driver directly reconstructs the open-square Hamiltonian and
diagonalizes it at every SCF iteration; it does not form its reference from an
MPO. It writes under the separate campaign name
`separable_aubry_andre_lside5_seed0p1_dense_hf`, so both calculation paths
retain their densities, bond orders, SCF histories, energy components, and
provenance independently.

## P1.2: fixed-Hamiltonian SP2 cap diagnostic

Before retrying the MPO SCF calculation, run the initial-SP2 diagnostic:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_2_fixed_sp2.slurm
```

This submits one job per input Hamiltonian. Each job runs the same first SCF
Hamiltonian, `H0 + S`, at `maxdim = 128, 256, 384, 512`, without updating any
Hartree or Fock field. Every run writes `iterations.csv` with the trace,
candidate polynomial traces, branch, rank, cap-contact flag, idempotency, and
per-step resources. `summary.toml` records the production SP2 stopping result.

The result root is
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside5_seed0p1_fixed_sp2/`.
Interpret this experiment before changing a field extractor: if SP2 fails in
this isolated problem after reaching the cap, the limitation is purification
compression rather than SCF feedback.

## P1.3: tight fixed-Hamiltonian interval control

The dense reference gives the clean initial spectrum as
`[-3.983143..., 3.983143...]`. Run this single control before altering any
production SCF bound:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_3_tight_interval_sp2.slurm
```

It changes only the fixed initial-SP2 interval from `[-12.5, 12.5]` to the
safe dense enclosure `[-4.25, 4.25]`, retaining `maxdim=256`,
`itensors_tol=1e-14`, and every other P1.2 setting. This interval is only
valid for `H0 + S`; it is deliberately not proposed as an SCF production
bound after Hartree and Fock fields have been added.

## P1.4 and P1.5: independent tight-scaling controls

After submitting P1.3, these array controls can run independently:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_4_tight_cap_ladder.slurm

sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_5_interval_ladder.slurm
```

P1.4 holds the interval at `[-4.25, 4.25]` and runs `maxdim=128,384`; P1.3
already supplies the `maxdim=256` middle point. P1.5 holds `maxdim=256` and
uses `[-4.05,4.05]`, `[-5,5]`, and `[-6,6]`; it omits the in-progress P1.3
`[-4.25,4.25]` case. All five jobs isolate one numerical control at a time.

## P1.6: branch-free McWeeny control

The clean fixed initial Hamiltonian is particle-hole symmetric at half filling,
so its exact chemical potential is `mu=0`. Run the direct-McWeeny comparison:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_6_mcweeny_mu_control.slurm
```

It uses the same `[-4.25,4.25]` linear initialization and `maxdim=256` as
P1.3, but replaces the SP2 branch recurrence with
`P -> 3P^2 - 2P^3`. Trace is reported as a symmetry invariant, not used to
select a polynomial or enforce the grand-canonical McWeeny occupation.

## P1.7: looser McWeeny compression control

To assess a practical `1e-3`-scale projector target, run the same P1.6
control with only the MPO compression cutoff changed from `1e-14` to `1e-8`:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_7_mcweeny_tol1e8.slurm
```

The diagnostic still records trace and idempotency at every step and retains
the production McWeeny `1e-6` stopping criterion. Thus a nonconverged result
can still demonstrate whether looser compression reaches a stable practical
residual near `1e-3` at lower cost.

## P1.8: fixed-\(\mu\) gap ladder

Run the three checkerboard gaps as independent array tasks:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_8_mcweeny_gap_ladder.slurm
```

The target gaps `0.5`, `1.0`, and `2.0` use checkerboard masses `0.25`,
`0.5`, and `1.0`, respectively. Every task holds `mu=0`, the interval
`[-4.25,4.25]`, `maxdim=256`, `itensors_tol=1e-8`, and 50 McWeeny steps
fixed. The output root is
`$MPO_RESULTS_ROOT/square_lside5_mcweeny_gap_ladder/`.

The ladder identifies gap `0.2` as the difficult case and gap `1.0` as a
transition case for the Horner-form comparison below.

## P1.9: standard versus Horner-form McWeeny

Submit the four controlled fixed-Hamiltonian jobs:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_9_mcweeny_horner_comparison.slurm
```

Tasks 1 and 2 compare the standard and Horner evaluation orders at gap `0.2`;
tasks 3 and 4 repeat the comparison at gap `1.0`. Every task uses `mu=0`,
the interval `[-4.25,4.25]`, `maxdim=256`, `itensors_tol=1e-8`, and 50
iterations. The two forms evaluate the same exact polynomial,
`3P^2-2P^3 = P^2(3I-2P)`, but MPO truncation occurs in a different order.

Each result writes full-update wall time, allocations, the intermediate bond
dimension, the next-density bond dimension, trace, and idempotency. Outputs
are stored under
`$MPO_RESULTS_ROOT/square_lside5_mcweeny_horner_comparison/`.

## P1.10: best-iterate McWeeny selection

Fixed-\(\mu\) McWeeny now retains the evaluated MPO with the lowest
idempotency residual. If the `1e-6` production threshold is not reached, the
solver returns that best iterate with `converged=false` instead of returning a
later deteriorated state. `PurificationResult.selected_iteration` and the
`purification_selected_iteration` SCF-history column identify the chosen
state; `iterations` still records the total work performed.

Horner remains opt-in for an SCF campaign:

```julia
mcweeny_form = :horner
```

A canonical trace guard is also opt-in and requires both fields:

```julia
mcweeny_trace_target = 512.0
mcweeny_trace_tolerance = 1e-3
```

The target is deliberately not inferred from `params.density`, because
unconstrained fixed-\(\mu\) McWeeny must remain grand-canonical. A
best-iterate result that did not meet `1e-6` remains nonconverged; a campaign
must explicitly set `allow_unconverged_purification=true` if it intends to
continue SCF with a separately justified practical tolerance.

## P1.11: large-gap SP2 control

Test canonical SP2 on a fixed clean Hamiltonian with checkerboard mass
`S(x,y)=2(-1)^(x+y)`. Its HOMO and LUMO lie near `-2` and `+2`, respectively,
so the total one-particle gap is `4`:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_11_sp2_gap4.slurm
```

The run uses `maxdim=256`, `itensors_tol=1e-14`, 50 SP2 steps, half filling,
and the safe fixed-Hamiltonian interval `[-4.75,4.75]`. It performs no SCF
and writes to `$MPO_RESULTS_ROOT/square_lside5_sp2_gap4/maxdim_256/`.

## P3.1: first CUDA SP2 calibration

Before moving a full SCF iteration to the GPU, repeat the validated fixed
gap-2 Hamiltonian on one accelerator:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p3_1_cuda_sp2_gap2.slurm
```

This uses the operational calibration (`maxdim=256`,
`itensors_tol=1e-10`, SP2 idempotency tolerance `2e-4`, and relative trace
tolerance `1e-6`). The Hamiltonian and density MPO are transferred to CUDA;
the complete SP2 recurrence, including truncating factorizations, runs there.
Only the final projector returns to the CPU for comparison with exact
diagonalization.

The submission loads MN5's NVIDIA HPC SDK, reads its CUDA version from
`nvcc`, and provisions CUDA.jl's `local_toolkit` preference before starting a
fresh Julia process. This is required because the transferred offline depot
was prepared without an NVIDIA driver and therefore cannot select a CUDA
runtime during precompilation.

The result directory is
`$MPO_RESULTS_ROOT/square_lside5_sp2_cuda_gap2/maxdim_256_tol_1e-10/`.
Inspect `summary.toml`, `sp2_progress.txt`, and `cuda_status.txt`. Do not use
the CUDA path for SCF until convergence, density and bond errors, energy
error, and iteration count agree with the established CPU calibration.

## P1.12: SP2 gaps 1 and 2 with cap control

Run gaps `1` and `2` at both `maxdim=256` and `384`:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p1_12_sp2_gap1_gap2_cap_ladder.slurm
```

The four array tasks hold the cutoff at `1e-14`, use 50 SP2 steps, half
filling, and the common safe interval `[-4.25,4.25]`. The output hierarchy is
`$MPO_RESULTS_ROOT/square_lside5_sp2_gap1_gap2/<gap>/<maxdim>/`. Comparing
the two caps distinguishes a gap-limited projector from a cap-limited SP2
trajectory.

## P2.1: gap-1 SP2 versus exact diagonalization

Compare the converged gap-1 SP2 projectors at both caps with an independently
constructed dense `1024 x 1024` occupied projector:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p2_1_sp2_dense_gap1.slurm
```

The two array tasks use `maxdim=256` and `384`. Each reconstructs and
diagonalizes the fixed Hamiltonian directly, then records all site-density
and nearest-neighbor bond-order errors, checkerboard order, one-body energy,
trace, idempotency, Hermiticity, and stationarity. The MPO is not converted
to a full dense matrix. Outputs are written below
`$MPO_RESULTS_ROOT/square_lside5_sp2_dense_gap1/`.

## P2.2: strict gap-1 SP2 stopping threshold

Test whether the `maxdim=384` gap-1 projector continues improving when the
SP2 idempotency threshold is tightened from the default `1e-3` to `1e-4`:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p2_2_sp2_dense_gap1_strict.slurm
```

This repeats the complete exact-diagonalization comparison and writes to
`$MPO_RESULTS_ROOT/square_lside5_sp2_dense_gap1_strict/`. It is intentionally
a single `maxdim=384` run: the question is whether additional SP2 iterations
improve the projector before fixed-cap truncation produces stagnation.

## P2.3: exact and SP2 projector compression

For the checkerboard potential `S=±1` (gap 2), compare the intrinsic QTT
compression of the exact dense occupied projector with post-compression of
SP2 projectors generated at `maxdim=256` and `384`:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p2_3_projector_compression_gap2.slurm
```

The three parallel tasks are exact compression with `maxdim=1024`, SP2 with
`maxdim=256`, and SP2 with `maxdim=384`. Each uses the cutoff ladder
`1e-14,1e-12,1e-10,1e-8,1e-6`. Every candidate is converted back to a dense
matrix and checked against exact diagonalization for global projector error,
trace, idempotency, Hermiticity, stationarity, eigenvalue bounds, all site
densities, nearest-neighbor bonds, and one-body energy. Results are written
under `$MPO_RESULTS_ROOT/square_lside5_projector_compression_gap2/`.

## P2.4: SP2 with an operational `1e-8` cutoff

Test whether SP2 can follow the intrinsically low-rank gap-2 projector
directly when every MPO operation uses `itensors_tol=1e-8`, rather than
applying that cutoff only after convergence:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p2_4_sp2_cutoff_1e8_gap2.slurm
```

The two array tasks use `maxdim=128` and `256`. Both run trace-correcting SP2
on the same `S=±1` gap-2 Hamiltonian, compare the resulting projector with
exact diagonalization, and record a final `1e-8` recompression control.
Outputs are written below
`$MPO_RESULTS_ROOT/square_lside5_sp2_operational_tol_1e8_gap2/`.

## P2.5: proposed SP2 production error budget

Calibrate the proposed hierarchy below the eventual `1e-3` SCF target on the
fixed gap-2 Hamiltonian:

```text
operational itensors_tol       = 1e-10
SP2 idempotency tolerance      = 2e-4
SP2 relative trace tolerance   = 1e-6
final recompression cutoff     = 1e-8
```

Submit the `maxdim=128` and `256` tasks with:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_p2_5_sp2_error_budget_gap2.slurm
```

The relative trace tolerance is converted to the size-aware absolute
tolerance `Ne*1e-6`; it is used only for convergence acceptance. SP2 branch
selection retains its original numerical tolerance. Results are written to
`$MPO_RESULTS_ROOT/square_lside5_sp2_error_budget_gap2/`.

## Zip-up versus variational square compression

`submit_sp2_variational_square.slurm` is an isolated fixed-Hamiltonian
diagnostic for the quasiperiodic `V2=0.5` case. It does not run SCF. At each
SP2 step it builds a `maxdim=1024` reference square, then compares:

- the existing one-pass `apply(rho, rho)` zip-up truncation;
- a two-sweep global variational compression initialized from that zip-up
  square and fitted to the explicit high-cap reference;
- a two-sweep implicit variational compression of `rho*rho`, represented as
  an MPO acting on a vectorized MPO without constructing the high-cap product.

The two array tasks use target caps 256 and 512. Read `iterations.csv` to
compare `square_relative_error` and `relative_exact_projector_error` between
the two independently advanced trajectories. This experiment is deliberately
restricted to `L_side=5`; the high-cap reference square is a validation tool,
not a proposed production algorithm.

Before attempting that full comparison on an accelerator, validate the
variational fitting paths with three small GPU iterations:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_sp2_variational_square_cuda_smoke.slurm
```

This requests one H100 and uses target `maxdim=256`, reference `maxdim=512`,
two fitting sweeps, and three SP2 iterations. It records both a one-second
`nvidia-smi` memory trace and synchronized CUDA timings. Passing this smoke
test establishes backend coverage only; it does not establish convergence or
justify replacing production zip-up multiplication.

After the Float64 smoke test reproduces the CPU branch sequence and
compression errors, submit the complete two-cap GPU comparison with:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_sp2_variational_square_cuda.slurm
```

This runs 20 fixed-H iterations at target caps 256 and 512, each against a
reference cap of 1024. Results and one-second GPU-memory traces are placed
under
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside5_sp2_square_fit_cuda_float64/`.

## Implicit-only fixed-H SP2

`submit_sp2_implicit_square_cuda.slurm` runs the `v2_0p5_u_1` fixed initial
Hamiltonian on one H100 for target bond dimensions 256 and 512. It advances a
single SP2 trajectory for at most 20 iterations using a two-sweep implicit
variational fit for every square. It does not construct an exact projector or
a higher-cap reference square, and it uses the current density MPO as the
initial fit so that no explicit MPO–MPO square is required.

The operational tolerances are cutoff `1e-8`, relative trace error `1e-6`, and
idempotency residual `2e-4`. Per-iteration fit/update timings, bond dimensions,
and free device memory are written to `iterations.csv`; Slurm also samples
whole-process GPU memory once per second.

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_sp2_implicit_square_cuda.slurm
```

## Sparse-vector KPM local-observable diagnostic

`submit_kpm_local_cuda.slurm` tests a different representation: it never
constructs an MPO density matrix. For the fixed initial `v2_0p5_u_1`
Hamiltonian, it applies a Jackson-damped Chebyshev approximation to blocks of
probing vectors and estimates only `rho[i,i]` and open-boundary
nearest-neighbour `rho[i,j]`. Exact diagonalization supplies the independent
reference. This is not an SCF run.

Array tasks 1--3 use all 1024 Hadamard probes at 128, 256, and 512 moments.
Because `Z*Z'/1024=I`, these runs have no probing error and isolate the KPM
polynomial error. Tasks 4--9 hold 512 moments and compare randomized Hadamard
and Rademacher sets with 32, 128, and 512 probes. They measure the extra local
probing error relevant to a scalable calculation.

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_kpm_local_cuda.slurm
```

Each task writes `summary.toml`, `density.csv`, `bonds.csv`, a process resource
report, and a one-second GPU-memory trace under
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside5_kpm_local/`. The supplied
`[-4.1,4.1]` interval is validated against the exact fixed-H spectrum before
the recurrence starts. Density and bond tables include the across-probe sample
variance, variance divided by the absolute estimated mean, coefficient of
variation, and naive relative standard error. The latter is a statistical
standard error for independent Rademacher probes; Hadamard columns are
correlated, and the complete Hadamard reconstruction is instead exact by
orthogonality. `probe_statistics.csv` also records the error distribution of
every individual probe against exact diagonalization.

## Complete KPM-SCF versus dense HF

The strict comparison runs the complete self-consistent `V2=0.5`, `U=1`
calculation twice on the same `32 x 32` open square: once with KPM local
observables and once with independent dense diagonalization. The KPM path uses
all `N=1024` Hadamard vectors, so this comparison has polynomial error but no
probing error. Both solvers require two consecutive iterations below the
relative SCF tolerance `2e-4`.

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside5/submit_kpm_dense_strict_comparison.slurm
```

The job writes both complete result trees and a key-aligned
`comparison.toml` under
`$MPO_RESULTS_ROOT/kpm_dense_lside5_strict/job_<job-id>/`.
