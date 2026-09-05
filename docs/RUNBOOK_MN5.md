# MN5 runbook

## Environment

Submit from the repository root and keep results external:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
export JULIA_DEPOT_PATH=/gpfs/projects/epor78/cluster_depot
```

Launchers load Julia 1.12.6. CUDA launchers load the MN5 NVIDIA HPC SDK, derive `CUDA_PATH` from `nvcc`, and initialize CUDA.jl's local-runtime preference before running the solver. This restart is deliberate: CUDA.jl may have been precompiled on a login node without a driver. Test CUDA only on an allocated GPU node. The supplied H100 launchers request one GPU and 20 CPUs.

## Submission and monitoring

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside6/submit_public_mpo_cuda.slurm

squeue -j "$JOB_ID" -o "%.18i %.10T %.12M %.6C %R"
sacct -j "$JOB_ID" --format=JobID,State,ExitCode,Elapsed,MaxRSS
```

Use the result path printed by the job, not a guessed `job_JOBID` path. Public campaigns write `input.toml`, `metadata.toml`, `observables.toml`, `scf_history.csv`, `site_density.csv`, and `bond_order.csv`. KPM expert runs also write `final_audit.toml`.

## Synchronization and provenance

Inspect the sync file list before using deletion. Never put generated results inside the source repository. Retain campaign/input snapshots, Slurm stdout and stderr, resource records, full result directories, and the local source revision. The cluster may lack Git, so capture the revision locally before sync or preserve the copied campaign source with each result.
