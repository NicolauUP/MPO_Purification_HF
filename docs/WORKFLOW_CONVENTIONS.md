# Workflow conventions and migration boundary

This document prevents old experiment launchers from becoming templates for
new work. It classifies the execution and output conventions that currently
coexist in the repository; it does not invalidate the scientific provenance
of historical campaigns.

## Canonical entry points

Use these for new production calculations:

- `bin/mpohf.jl preflight|run` for dense, MPO, and KPM public campaigns;
- `CampaignSpec` containing explicit `CaseSpec` values for physical and
  numerical configuration;
- `runs/submit_public_campaign.slurm` as the cluster-neutral CPU launcher;
- `runs/template/submit_mpo_cuda.slurm` as the MareNostrum5 MPO/CUDA launcher;
- `runs/common/slurm_cuda_runtime.sh` for the shared MN5 CUDA runtime setup;
- `bin/mpohf_postprocess.jl` for observables deferred from a converged MPO run.

The canonical result directory is immutable and has the form

```text
<output-root>/<campaign>/task_<index>_<label>/
```

Its standard files are `campaign.jl`, `input.toml`, `metadata.toml`,
`observables.toml`, `scf_history.csv`, `site_density.csv`, and
`bond_order.csv`. A converged MPO run also writes `converged_state.h5`.
The campaign source is the executable input; `input.toml` is a resolved,
descriptive summary because Julia functions cannot be reconstructed from TOML.

## Supported configuration schema

New campaign sources must bind `campaign` to a `CampaignSpec`. Keep four
concerns separate:

1. `ChainModel`, `SquareModel`, or `GraphModel`: physical model and geometry;
2. `QTTSettings`: representation and compression;
3. `SCFSettings` and `KPMSettings`: numerical algorithm controls;
4. `RuntimeSettings` or CLI options: CPU/CUDA execution choice.

Cluster resources and module initialization belong in SLURM templates, not in
campaign sources. Physics and convergence controls belong in campaign sources,
not in SLURM files.

## Historical and expert-only conventions

The following are preserved for reproducibility, but are not starting points
for new production work:

- tuple campaigns shaped as `(name=..., runs=[...])`;
- `runs/common/execute_1d_campaign.jl`,
  `execute_square_campaign.jl`, and the separate dense variants;
- launchers that configure `CUDA.set_runtime_version!` inline or source
  `slurm_cuda_job_depot.sh`;
- fixed-field SP2, KPM estimator, replay, compression, profiling, and audit
  scripts in `runs/common/` or `benchmark/`;
- result trees organized around `job_$SLURM_JOB_ID`, `maxdim_*`, `side_*`, or
  script-specific filenames such as `density.csv` and `summary.toml`.

These files often encode a deliberately different method or diagnostic. Do not
bulk-convert them. Migrate one only when rerunning that exact experiment, after
checking its README, call path, parameter values, and expected artifacts.

## Checkpoint and restart semantics

`converged_state.h5` is currently a final-state checkpoint. It stores the MPO
Hamiltonian, fields, and density so that observables can be measured later.
It is written only after convergence and is not an SCF continuation or
mid-run recovery checkpoint. `read_result` is likewise an analysis/provenance
reader, not a restart API.

Files named snapshot, replay, or fixed-SP2 checkpoint belong to specific
diagnostics and must not be presented as general restart support. A future
restart feature needs an explicit format version, saved iteration/mixing
state, compatibility checks, and tests before it becomes canonical.

## Conservative migration recipe

For a legacy campaign that needs a new run:

1. Preserve the original campaign and launcher unchanged.
2. Translate one case at a time to `CampaignSpec`/`CaseSpec`, retaining every
   physical parameter, bound, tolerance, seed, method, and compression cap.
3. Run `bin/mpohf.jl preflight` and a deterministic small validation case.
4. Compare the translated configuration with the original before submission.
5. Use the canonical CPU or MN5 CUDA template and a new campaign name/label.
6. Record the legacy source path and comparison evidence in the campaign
   README. Never overwrite the historical result directory.

If a legacy script implements a capability absent from the public API, leave
it expert-only and document that gap instead of hiding it behind a nominally
canonical launcher.
