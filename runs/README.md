# Production campaigns

`runs/` contains tracked campaign definitions and submission templates. It
does not contain generated results. A campaign is a finite, explicit list of
parameter points; each Slurm array task selects one point and writes an
immutable result directory outside this repository.

## Workflow

1. Copy `runs/1d/template/` to a descriptively named campaign directory.
2. Replace the example `campaign.jl` entries with the parameter points to run.
3. Copy and edit `submit.slurm` for the target cluster and campaign size.
4. Set `MPO_RESULTS_ROOT` to an external storage location, then submit the
   array from the package root.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/<account>/MPO_HF_results
sbatch runs/1d/<campaign>/submit.slurm
```

The generic driver refuses to overwrite an existing result directory. Reruns
must use a new campaign name or an explicit new label, preserving provenance.

## Per-task output contract

For task `i`, the driver creates:

```text
$MPO_RESULTS_ROOT/<campaign>/task_<i>_<label>/
  input.jl             Snapshot of the exact campaign definition
  selection.toml       Selected task, bounds, method, and label
  metadata.toml        Git/Julia/Slurm provenance and final status
  scf_history.csv      One row per SCF iteration
  observables.toml     Final scalar observables and numerical diagnostics
  site_density.csv     Final local density
  bond_order.csv       Final nearest-neighbour bond order
```

The first version deliberately writes no MPO checkpoint. Checkpoints should
be added only after deciding the stable file format and restart semantics.

`input.jl` is copied rather than reconstructed from TOML because a run may
contain Julia functions for hopping, potential, or seeds. Together with the
selected array index and Git revision it is the executable scientific input.

## What the driver does not decide

It does not infer spectral bounds, choose a purification method, or decide
whether a result is physically acceptable. Those choices remain explicit in
each `campaign.jl`. A completed Slurm task may still report unconverged SCF;
inspect `metadata.toml`, `scf_history.csv`, and `observables.toml`.
