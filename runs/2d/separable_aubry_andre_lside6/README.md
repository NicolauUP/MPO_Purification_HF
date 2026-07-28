# 64x64 separable quasiperiodic pilot

This campaign fixes one physical case:

- `side_level = 6`, so the square is `64 x 64`, `N = 4096`;
- `V2 = 0.5`, `U = 1.0`, half filling;
- `tau = sqrt(2) - 5/6`;
- separable hopping modulation in `x` and `y`;
- checkerboard initial seed of amplitude `0.1`;
- operational MPO cutoff `1e-10`.

The fixed initial Hamiltonian uses `[-6.35, 6.35]`: its hopping row sum is at
most `6`, the seed contributes `0.1`, and the remaining `0.25` is margin. The
full SCF jobs use `[-8.5, 12.5]`. For repulsive `U=1`, this includes
`0 <= V_H <= 4`, a Fock row-sum bound of `2` from
`|rho_ij| <= 1/2`, the hopping row sum `6`, the seed, and margin.

The three experiments answer different questions:

1. `submit_exact_projector_compression.slurm` diagonalizes the fixed initial
   Hamiltonian once and compresses its exact occupied projector at
   `maxdim = 128, 256, 512, 1024`, with cutoff `1e-8`.
2. `submit_fixed_sp2_cap_ladder_cuda.slurm` compares CUDA SP2 directly with
   exact diagonalization at `maxdim = 256, 512, 1024`.
3. `submit_full_scf_cap_ladder_cuda.slurm` runs the full Hartree--Fock SCF
   calculation at the same three caps.
4. `submit_fixed_sp2_truncation_cuda.slurm` keeps `maxdim=512`, uses the
   exact-informed fixed-H interval `[-4.1, 4.1]`, and compares operational
   cutoffs `1e-8` and `1e-10`. Nonconverged SP2 runs return their best measured
   iterate rather than an unmeasured final polynomial update.

Inspect experiments 1 and 2 before submitting experiment 3. In particular,
the `maxdim=1024` SP2 work MPO may exceed one H100's available memory even
when the final exact projector compresses below that cap.

Submit from the repository root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks

sbatch runs/2d/separable_aubry_andre_lside6/submit_exact_projector_compression.slurm
sbatch runs/2d/separable_aubry_andre_lside6/submit_fixed_sp2_cap_ladder_cuda.slurm
```

After reviewing those results, submit selected SCF tasks. For example,
start with `maxdim=256`:

```bash
sbatch --array=1 runs/2d/separable_aubry_andre_lside6/submit_full_scf_cap_ladder_cuda.slurm
```

Array indices `1,2,3` correspond to `maxdim=256,512,1024`.
