# Run definitions

This directory holds **tracked definitions** of scientific calculations. It is
not a location for cluster outputs, copied Julia depots, checkpoints, or Slurm
logs. Keeping definitions separate from results makes it possible to identify
exactly which code revision and numerical settings generated a result.

## Layout

```text
runs/
  1d/                 Reproduction cases and parameter scans for chains
  2d/                 Square-lattice development and production definitions
  templates/          Starting-point notes for a new run definition
  README.md           This policy
```

When a calculation becomes important, give it a dedicated tracked directory,
for example:

```text
runs/1d/reference_modulation_a/
  definition.jl       Parameters, spectral bounds, and method choice
  submit.slurm        Resource request and invocation
  README.md           Physical purpose and expected invariants
```

The cluster result directory should instead live outside the repository, for
example:

```text
/gpfs/projects/epor78/MPO_HF_runs/1d_reference_modulation_a/<run-id>/
```

This prevents `rsync --delete` or a source-tree cleanup from destroying
results, and prevents large outputs from being committed accidentally.

## Minimum record for every scientific run

Before submission, record in `definition.jl` or the run README:

1. Physical model: geometry, `L`, hopping, `U`, `W`, filling, boundary
   condition, and initial seed `S`.
2. Numerical settings: `tci_tol`, `itensors_tol`, `itensors_maxdim`,
   purification method and steps, SCF mixing/tolerance/iteration limit, and
   explicit `(H_min, H_max)` spectral bounds.
3. Intended validation: e.g. a known 1D reference, non-interacting limit,
   dense comparison, or a symmetry/particle-number invariant.
4. The exact Git commit and Julia version used on the cluster.

The result directory should contain, at minimum, the copied definition,
Slurm stdout/stderr, a scalar SCF history, final observables, and metadata
describing the host, Julia version, depot, and Git revision. The package does
not yet provide one unified production-output writer; this is the target
contract for the first run driver rather than a claim about present output.

## Recommended workflow

1. Start with a small 1D definition under `runs/1d/` and reproduce a result
   using both MPO and a dense small-system reference where possible.
2. Freeze that definition once its density, bond order, energy, particle
   number, and stationarity are understood.
3. Reuse the same output contract for the first `4×4` square calculation.
4. Move to larger systems only after checking purification convergence,
   bond-dimension growth, and SCF history for the smaller case.

Do not use a result merely because a job exited successfully. Record and
inspect the purification status, `Tr(rho)`, idempotency, Hermiticity,
commutator residual, and final bond dimension.
