# 256x256 CDW seed-stability campaign

This nine-task GPU array tests whether the symmetry-broken Hartree--Fock
solution is selected independently of the magnitude and sign of an explicit
checkerboard seed. It does **not** repeat the exactly uniform seed, because an
exactly symmetry-preserving initial field can remain in that sector even when
the CDW solution has lower energy.

For each representative point it runs checkerboard seed amplitudes
`+0.1`, `-0.1`, and `+1e-3`:

| Label | \(V_2\) | \(U\) | production \(M\) | audit \(M\) |
|---|---:|---:|---:|---:|
| `v2_0_u_0p5` | 0.0 | 0.5 | 1200 | 1600 |
| `v2_0p5_u_0p5` | 0.5 | 0.5 | 1600 | 2000 |
| `v2_2_u_0p3` | 2.0 | 0.3 | 2400 | 3200 |

All jobs use the validated safeguarded Pulay policy (history 3, warm-up 3,
regularization \(10^{-12}\), coefficient cap 16, full Pulay damping and 0.5
linear fallback), up to 60 SCF iterations. The independent final audit uses
2048 Hadamard probes. `final_audit.toml` also records the audited checkerboard
order.

Submit from the repository root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside8/submit_kpm_seed_stability.slurm
```

Interpretation requires all three seeds to converge. For the clean bipartite
case, opposite checkerboard signs can be symmetry-related and hence should
have equal energy and opposite order. The weak seed tests whether the CDW
attractor is reached without imposing a large initial field.
