# 256x256 KPM moment-convergence campaign

This GPU array selects a common production KPM resolution before a broad
interacting \(U,V_2\) phase map. It runs three representative points:

| Label | \(V_2\) | \(U\) |
|---|---:|---:|
| `v2_0_u_0p5` | 0.0 | 0.5 |
| `v2_0p5_u_0p5` | 0.5 | 0.5 |
| `v2_2_u_0p3` | 2.0 | 0.3 |

Each point is evaluated from checkerboard and uniform initial fields at
\(M=400,800,1200,1600\), with \(R=1024\) production Hadamard probes. Each
completed SCF is independently audited at \(M=2000\), \(R=2048\).

The campaign uses the validated safeguarded Pulay policy: three-vector
history, three-step warm-up, regularization \(10^{-12}\), coefficient cap 16,
unit Pulay damping, and 0.5 linear fallback mixing.

Submit from the repository root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside8/submit_kpm_moment_convergence.slurm
```

The array limits concurrent tasks to six GPUs. Results are grouped by physical
point, seed, and production moment count beneath one array-job directory.
Compare both seeds and require agreement of production observables with the
independent audit before choosing a common \(M\) for a phase map.
