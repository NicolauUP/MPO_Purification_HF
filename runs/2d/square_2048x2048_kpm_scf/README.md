# KPM-HF at 4,194,304 sites

This is the first one-H100 size extrapolation after the byte-for-byte
reproducible `1024 x 512` KPM-SCF validation. It runs a `2048 x 2048` open
square with the same half-filled separable quasiperiodic hopping, `V2=0.5`,
`U=1`, checkerboard seed, KPM degree `M=1200`, and 1024 coordinate-interleaved
Hadamard probes. The independent final audit uses `M=1600` and 2048 probes.

The calculation uses the allocation-free Float64 CUDA Chebyshev recurrence.
The 524k validation showed that it gives byte-identical physical output while
reducing peak GPU memory from 64,950 MiB to 5,302 MiB.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/square_2048x2048_kpm_scf/submit_full_kpm_scf_cuda.slurm
```

Expected runtime is roughly 4--5 hours if SCF converges in about 20 iterations.
The output directory is
`$MPO_RESULTS_ROOT/square_2048x2048_kpm_scf/job_<job-id>/`.

Accept the result only after checking `scf_converged`, the final independent
audit differences, trace error relative to `N_e`, and the sampled GPU-memory
peak.
