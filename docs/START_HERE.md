# Start here

This package studies zero-temperature, spinless nearest-neighbour Hartree--Fock models with three complementary solvers:

| Method | Best use | Current scope |
|---|---|---|
| `:dense` | exact small-system reference | open 1D and equal open squares |
| `:mpo` | QTT/MPO purification research | equal open squares; CUDA path is hybrid |
| `:kpm` | large GPU calculations | open power-of-two squares, rectangles, and sparse graphs |

Choose the method from the scientific question, not apparent speed. Every new model should first be checked with dense diagonalization where feasible.

## Recommended public workflow

Work from the repository root. A public campaign is Julia source that binds a `CampaignSpec`; it owns the physics and numerical controls. The CLI selects a case and an explicit backend.

```bash
julia --startup-file=no --project=. bin/mpohf.jl preflight \
  runs/2d/separable_aubry_andre_lside6/campaign_public.jl \
  --task 1 --method mpo --backend cuda

julia --startup-file=no --project=. bin/mpohf.jl run \
  runs/2d/separable_aubry_andre_lside6/campaign_public.jl \
  --task 1 --method mpo --backend cuda \
  --output-root "$MPO_RESULTS_ROOT" --verbose all
```

Run `preflight` first. It rejects unsupported method/model combinations and reports the actual runtime path. MPO CUDA purification is currently hybrid: Hartree/Fock extraction remains host-side. CUDA KPM is more device-oriented.

## Read in this order

1. [Scientific conventions](SCIENTIFIC_CONVENTIONS.md).
2. [MN5 runbook](RUNBOOK_MN5.md).
3. [Validation matrix](VALIDATION_MATRIX.md).
4. [Current status](CURRENT_STATUS.md).
5. [Public API migration](plans/public_api_redesign.md).
6. [Workflow conventions and migration boundary](WORKFLOW_CONVENTIONS.md).

Campaign READMEs are the authority for campaign-specific provenance. Legacy `runs/common/` and `benchmark/` scripts are expert diagnostics, not defaults.

## Minimum acceptance checklist

Check `scf_converged`, termination reason, trace, idempotency, SCF residuals, two-cycle status, seed sensitivity, and the applicable dense comparison or independent KPM audit. A completed Slurm job is not itself a physical result. Results are immutable: create a new label for a rerun rather than overwriting numerical provenance.
