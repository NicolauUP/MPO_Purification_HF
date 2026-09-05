# Selected-entry KPM scaling at 256 x 256

This fixed-H diagnostic uses `L=16`, so `L_side=8`, the open square is
`256 x 256`, and `N=65,536`. It retains the `V2=0.5`, `U=1`, checkerboard
seed `0.1` case, but performs neither SCF nor exact diagonalization.

The safe interval is `[-6.35,6.35]`. The run uses 800 Jackson--Chebyshev
moments, approximately preserving the spectral resolution of 512 moments on
`[-4.1,4.1]`.

Partial Hadamard estimates at `R=256,512,1024,2048,4096` are compared with
deterministic entries of the same KPM polynomial. Those entries are generated
by applying KPM to basis vectors for an `8 x 8` coordinate grid and its valid
right/up bonds. Counts 512, 1024, and 2048 are repeated with another randomized
sign seed. Therefore reported selected-entry differences are probing errors,
not polynomial errors.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside8/submit_kpm_hadamard_selected_cuda.slurm
```

Outputs are stored below
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside8_kpm_selected/`.

## All-site nested-probe convergence

The all-site run evaluates all 65,536 densities and all 130,560 open-boundary
nearest-neighbour bonds. It uses coordinate-aware nested Hadamard prefixes
with `R=512,1024,2048,4096,8192`; the last is the operational reference for
probing convergence at the fixed 800-moment KPM polynomial. Independent seeds
are also run at `R=1024,2048,4096`.

This is not an ED validation and does not measure the Jackson--Chebyshev
polynomial error. It tests whether the required probe count remains stable
when increasing the system from 8,192 to 65,536 sites. The output includes
all-site density, Hartree-field, horizontal-Fock, and vertical-Fock
differences, together with trace, energy, timing, and GPU-memory records.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside8/submit_kpm_hadamard_all_nested_cuda.slurm
```

Results are written below
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside8_kpm_all_nested/job_JOBID/`.

## Full KPM Hartree--Fock pilot

The full-SCF pilot closes the nearest-neighbour Hartree--Fock loop using only
KPM estimates of every density and real horizontal/vertical bond order. It
uses fixed coordinate-aware Hadamard probes throughout SCF, so the iterative
map is deterministic.

Production iterations use `M=1200`, `R=1024`. A trace-moment pass solves for
the half-filled chemical potential before each local-observable pass.
Per-iteration Gershgorin bounds rigorously enclose the current sparse
effective Hamiltonian. Convergence requires two consecutive iterations below
the existing `0.1%` relative field/density/bond threshold.

After SCF stops, an independent fixed-field audit uses `M=1600`, `R=4096` and
a different Hadamard sign seed.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside8/submit_full_kpm_scf_cuda.slurm
```
