# Production campaigns

`runs/` contains tracked campaign definitions and submission templates. It
does not contain generated results. A campaign is a finite, explicit list of
parameter points; each Slurm array task selects one point and writes an
immutable result directory outside this repository.

## Canonical workflow

1. Copy `runs/template/campaign.jl` to a descriptively named campaign
   directory and keep it in the public `CampaignSpec`/`CaseSpec` schema.
2. Replace the example case with explicit parameter points. Do not generate an
   implicit grid: every submitted point should be visible in version control.
3. Run `bin/mpohf.jl preflight` for the selected task, method, and backend.
4. Use `runs/template/submit_cpu.slurm` for CPU work or
   `runs/template/submit_mpo_cuda.slurm` for an MN5 MPO/CUDA run. Change
   cluster resources in a thin profile; do not move physics into SLURM.
5. Set `MPO_RESULTS_ROOT` to an external storage location and submit from the
   package root.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/<account>/MPO_HF_results
sbatch --export=ALL,\
MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT",\
MPOHF_CAMPAIGN=runs/1d/<campaign>/campaign.jl,\
MPOHF_METHOD=mpo \
  runs/template/submit_cpu.slurm
```

The generic driver refuses to overwrite an existing result directory. Reruns
must use a new campaign name or an explicit new label, preserving provenance.

## Per-task output contract

For task `i`, the driver creates:

```text
$MPO_RESULTS_ROOT/<campaign>/task_<i>_<label>/
  campaign.jl          Snapshot of the exact executable campaign definition
  input.toml           Resolved descriptive configuration
  metadata.toml        Git/Julia/Slurm provenance and final status
  scf_history.csv      One row per SCF iteration
  observables.toml     Final scalar observables and numerical diagnostics
  site_density.csv     Final local density
  bond_order.csv       Final nearest-neighbour bond order
```

Converged MPO runs also write `converged_state.h5`. This is a final-state
artifact for deferred measurement, not a mid-run restart mechanism.

`campaign.jl` is copied because a run may contain Julia functions for hopping,
potential, or seeds. `input.toml` is a self-describing resolved summary, not an
executable reconstruction. Together with the selected task and source revision
they provide the run provenance.

See [`../docs/WORKFLOW_CONVENTIONS.md`](../docs/WORKFLOW_CONVENTIONS.md) for
the canonical-versus-legacy boundary. Historical tuple campaigns, specialized
diagnostics, and their result layouts remain supported for provenance but are
not templates for new production work.

## What the driver does not decide

It does not infer spectral bounds, choose a purification method, or decide
whether a result is physically acceptable. Those choices remain explicit in
each `campaign.jl`. A completed Slurm task may still report unconverged SCF;
inspect `metadata.toml`, `scf_history.csv`, and `observables.toml`.
