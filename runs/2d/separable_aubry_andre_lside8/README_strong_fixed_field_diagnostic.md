# Strong-quasiperiodicity fixed-field KPM diagnostic

This one-GPU diagnostic separates KPM polynomial degree from Hadamard-probe
count for the (V_2=2, U=0.3), 256x256 CDW solution. It performs no SCF.
Instead, it reconstructs an effective Hamiltonian from the saved **audited**
density and bond order of a converged `checkerboard_plus` production run.

It evaluates four nested-Hadamard configurations:

| Configuration | Purpose |
|---|---|
| `M4000_R2048`, `M4000_R4096` | calibrated baseline and probe-count check |
| `M6000_R2048`, `M6000_R4096` | first degree refinement and probe-count check |
| `M8000_R2048`, `M8000_R4096` | conservative high-degree candidate and probe-count check |

The comparison is faithful to the saved final local observables, but the
historical SCF runner did not persist its mixed Hartree/Fock vectors. Thus the
field is reconstructed from audited observables rather than being bytewise the
same Hamiltonian used in the original last KPM step.

Submit after defining the converged source directory:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
export KPM_FIXED_FIELD_SOURCE="$MPO_RESULTS_ROOT/separable_aubry_andre_lside8_seed_stability/job_JOBID/v2_2_u_0p3/checkerboard_plus/M2400"

sbatch --export=ALL,MPO_RESULTS_ROOT,KPM_FIXED_FIELD_SOURCE \
  runs/2d/separable_aubry_andre_lside8/submit_kpm_strong_fixed_field_diagnostic.slurm
```

Read `summary.csv` for each configuration relative to the saved source and
`comparisons.csv` for the direct probe-order and polynomial-order differences.
