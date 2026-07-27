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
