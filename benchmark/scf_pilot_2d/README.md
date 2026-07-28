# Bounded 2D SCF pilot

This runner is for one scientifically inspectable square-lattice SCF case. It
is not a parameter sweep and intentionally limits `L_side` to 64 (`L_total=12`,
`N=4096`). Its purpose is to establish a stable SCF configuration before a
production campaign.

Each completed run writes:

- `metadata.toml` — full physical and numerical input, conservative spectral bounds, and project fingerprint;
- `progress.txt` — purification and SCF progress without terminal control codes;
- `scf_history.csv` — one row per SCF iteration, including residuals, energy, and MPO bond dimensions;
- `phase_timings.csv` — synchronized per-iteration timings for initialization,
  purification, density transfer to the host, CPU mean-field extraction, field
  transfer to the device, and device-side diagnostics;
- `observables.toml` — final particle number, checkerboard order, energy components, residuals, ranks, and five density/Hartree probes;
- `process_time.txt` — wall time, CPU time, allocations reported by `/usr/bin/time`.

Submit a uniform-seed pilot first:

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
  benchmark/scf_pilot_2d/scf_pilot_2d.slurm \
  --side-level 6 --U 0.3 --potential checkerboard --potential-amplitude 0.6 \
  --seed uniform --maxdim 256 --itensors-tol 1e-14 --steps 50 \
  --scf-mixing 0.5 --scf-tol 0.1 --scf-max-iterations 30 --padding 0.5
```

Repeat the identical case with `--seed checkerboard_plus` and
`--seed checkerboard_minus`. Compare the final energy, particle-number error,
stationarity residual, density order, and final ranks before interpreting any
seed dependence as a physical broken-symmetry state.

`--cpus-per-task=16` requests about 32 GiB under the current GPP policy while
keeping Julia single-threaded. This is memory headroom, not parallel speedup.
For a larger MPO cap, override it at submission, for example
`sbatch --cpus-per-task=32 ...` for roughly 64 GiB.

## H100 pilot

The CUDA pilot keeps SP2 purification and the effective-Hamiltonian diagnostics
on the GPU. Density matrices are transferred to the CPU for the current
Hartree/Fock extraction kernels, and the resulting fields are transferred back
to the GPU. The phase log measures this hybrid workflow directly.

Submit the calibrated `L_side=6` case with:

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
  benchmark/scf_pilot_2d/scf_pilot_2d_cuda.slurm
```

In addition to the common outputs, a successful CUDA run writes
`cuda_status.txt` with the selected runtime, device, and memory-pool state.
The submission configures CUDA.jl from the CUDA toolkit loaded on the allocated
compute node; it does not require internet access.

## Fixed-control size scaling

Four independent submissions extend the validated case without changing
\(\chi_{\max}=256\), truncation tolerances, physical parameters, seed, or SCF
controls:

| Script | `L_side` | `L_total` | Lattice sites |
|---|---:|---:|---:|
| `scf_scaling_lside7_cuda.slurm` | 7 | 14 | 16,384 |
| `scf_scaling_lside8_cuda.slurm` | 8 | 16 | 65,536 |
| `scf_scaling_lside9_cuda.slurm` | 9 | 18 | 262,144 |
| `scf_scaling_lside10_cuda.slurm` | 10 | 20 | 1,048,576 |

Submit them as separate jobs:

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
for level in 7 8 9 10; do
  sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
    "benchmark/scf_pilot_2d/scf_scaling_lside${level}_cuda.slurm"
done
```

Each job writes to
`$MPO_BENCHMARK_ROOT/scf_scaling_lsideLEVEL_cuda_JOBID`. Because the jobs are
independent, a timeout or memory failure at a larger size does not discard the
smaller-size measurements.
