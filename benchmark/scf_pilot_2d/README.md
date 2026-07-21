# Bounded 2D SCF pilot

This runner is for one scientifically inspectable square-lattice SCF case. It
is not a parameter sweep and intentionally limits `L_side` to 64 (`L_total=12`,
`N=4096`). Its purpose is to establish a stable SCF configuration before a
production campaign.

Each completed run writes:

- `metadata.toml` — full physical and numerical input, conservative spectral bounds, and project fingerprint;
- `progress.txt` — purification and SCF progress without terminal control codes;
- `scf_history.csv` — one row per SCF iteration, including residuals and MPO bond dimensions; the energy column is intentionally empty because the direct square-energy diagnostic is evaluated only after convergence;
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
