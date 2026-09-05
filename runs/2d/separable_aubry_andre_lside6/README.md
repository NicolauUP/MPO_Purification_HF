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
5. `submit_kpm_hadamard_ladder_cuda.slurm` keeps 512 Jackson--Chebyshev
   moments and tests sparse-vector local-observable extraction with
   `R=128,256,512,1024,2048,4096` Hadamard probes. The complete `R=N=4096`
   task has exactly zero probing error and isolates the polynomial error.
   Probe counts 512, 1024, and 2048 are repeated with a second randomized
   Hadamard sign pattern.

Inspect experiments 1 and 2 before submitting experiment 3. In particular,
the `maxdim=1024` SP2 work MPO may exceed one H100's available memory even
when the final exact projector compresses below that cap.

## Public production SCF entry point

The legacy diagnostic runners above remain useful for validation and profiling.
For the reproducible production (64\times64) MPO calculation, use the public
campaign instead:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
sbatch runs/2d/separable_aubry_andre_lside6/submit_public_mpo_cuda.slurm
```

The physical model, QTT settings, SP2 tolerances, Fock backend, and stopping
policy are all declared in `campaign_public.jl`. The Slurm wrapper only sets
the MN5 CUDA runtime and allocation. This configuration uses
`maxdim=512`, `cutoff=1e-10`, binary-carry Fock extraction, SP2 idempotency
tolerance `2e-4`, and relative SP2 trace tolerance `1e-6`. It intentionally
does not record the expensive *per-iteration* full energy on GPU; final
observables still include the converged total energy.

The result is written once to:

```text
$MPO_RESULTS_ROOT/separable_aubry_andre_lside6_public_mpo/task_0001_v2_0p5_u_1_chi_512
```

The writer refuses to overwrite this directory, preserving exact input and
runtime provenance.

Before changing the production case to DIIS/Pulay, run the lower-cap pilot:

```bash
MPOHF_CAMPAIGN=runs/2d/separable_aubry_andre_lside6/campaign_public_pulay_chi256.jl \
  sbatch runs/2d/separable_aubry_andre_lside6/submit_public_mpo_cuda.slurm
```

It keeps the physical model fixed and changes only the mixer to a guarded
four-entry Pulay history after three updates. Compare its SCF iteration count,
final residuals, final observables, and peak GPU memory with the established
linear χ=256 reference before using Pulay at χ=512.

Submit from the repository root:

```bash
export MPO_RESULTS_ROOT=/gpfs/projects/epor78/MPO_HF_results
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks

sbatch runs/2d/separable_aubry_andre_lside6/submit_exact_projector_compression.slurm
sbatch runs/2d/separable_aubry_andre_lside6/submit_fixed_sp2_cap_ladder_cuda.slurm
sbatch --export=ALL,MPO_RESULTS_ROOT="$MPO_RESULTS_ROOT" \
  runs/2d/separable_aubry_andre_lside6/submit_kpm_hadamard_ladder_cuda.slurm
```

After reviewing those results, submit selected SCF tasks. For example,
start with `maxdim=256`:

```bash
sbatch --array=1 runs/2d/separable_aubry_andre_lside6/submit_full_scf_cap_ladder_cuda.slurm
```

Array indices `1,2,3` correspond to `maxdim=256,512,1024`.
