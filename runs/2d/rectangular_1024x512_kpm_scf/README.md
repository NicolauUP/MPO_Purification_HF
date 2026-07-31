# Blocked KPM-HF at 524,288 sites

This production-scale pilot runs the full nearest-neighbour Hartree--Fock SCF
loop on a `1024 x 512` open rectangle (`N=524,288`) with the validated
`V2=0.5`, `U=1`, half-filled separable quasiperiodic model.

Production uses `M=1200`, `R=1024` in blocks of 128 probes. The final
independent fixed-field audit uses `M=1600`, `R=2048` with another seed.
Physical storage remains x-contiguous while Hadamard colors use interleaved
coordinate bits.

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/rectangular_1024x512_kpm_scf/submit_full_kpm_scf_cuda.slurm
```

This is a scaling experiment, not an independently validated solution. Accept
the result only if the SCF converges and the independent `M=1600`, `R=2048`
audit agrees with the production density and bond order at the chosen
tolerance. Outputs are written to
`$MPO_RESULTS_ROOT/rectangular_1024x512_kpm_scf/job_<job-id>/`; the adjacent
`_process_time.txt` and `_gpu_memory.csv` files record whole-job resources.

The allocation-free Chebyshev recurrence must first pass its CUDA/CPU
equivalence smoke test:

```bash
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/rectangular_1024x512_kpm_scf/submit_inplace_cuda_smoke.slurm
```

This verifies the CUDA `mul!` recurrence against the CPU implementation before
using it in another production-sized calculation.
