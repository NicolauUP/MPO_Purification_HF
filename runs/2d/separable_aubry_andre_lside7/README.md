# KPM validation against ED at 128 x 128

This campaign validates the fixed-H KPM local-observable estimator on a
`128 x 128` open square (`N=16,384`) with separable quasiperiodic hopping
`V2=0.5` and checkerboard seed `0.1`.

One exact diagonalization supplies the occupied-projector reference. The
eigenpairs also provide deterministic local observables for each finite KPM
polynomial, without constructing a dense projector. This separates:

- Jackson--Chebyshev polynomial error for `M=400,800,1200,1600`;
- Hadamard probing error for `R=512,1024,2048,4096`;
- total KPM error relative to ED.

Counts `R=1024,2048,4096` are repeated with an independent sign seed. Every
site density and every horizontal and vertical nearest-neighbour bond enters
the reported error norms.

Submit from the `MPO_Purification_HF` repository root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside7/submit_kpm_ed_decomposition_cuda.slurm
```

The main outputs are `polynomial_errors.csv`, `probing_errors.csv`, and
`metadata.toml` below
`$MPO_RESULTS_ROOT/separable_aubry_andre_lside7_kpm_ed/job_JOBID/`.
