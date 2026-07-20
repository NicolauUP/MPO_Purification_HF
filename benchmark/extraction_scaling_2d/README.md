# Square-lattice extraction benchmark

This benchmark measures the implemented open-square mean-field extractors on
`2^L_side × 2^L_side` lattices. Thus `L_side=10` means `L_total=20` QTT bits
and `N=1,048,576` sites.

It is deliberately not a full SCF calculation. Density preparation happens
outside the timed region. The timed kernels are:

- `hartree_tci` — current TCI-based square Hartree field;
- `hartree_tensorial` — experimental translation-based Hartree field;
- `fock_horizontal` and `fock_vertical` — the two implemented TCI Fock
  components, measured separately.

There is currently no square binary-carry Fock implementation, so this is not
a comparison of that 1D prototype. The two Hartree implementations are compared
as MPOs; Fock has no independent large-system implementation yet and is
reported for runtime, allocations, Hermiticity, and sampled physical bonds.

## Cluster preflight

Submit the preflight first. It runs only smooth synthetic densities up to a
`4×4` lattice and makes no large allocation.

```bash
cd /gpfs/projects/epor78/MPO_HartreeFock_Purification/MPO_Purification_HF
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
mkdir -p "$MPO_BENCHMARK_ROOT"

sbatch benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --side-levels 2 --sources smooth --warmups 1 \
  --repetitions-small 1 --repetitions-large 1
```

## Full scaling run

The default grid is `L_side=2,3,...,10`, equivalent to total QTT bit counts
`L_total=4,6,...,20`. The smooth source is always suitable for an extraction
scaling study. `sp2_gapped` additionally includes SP2 density construction;
it is useful only after a small preflight establishes that the chosen
truncation settings are feasible on the target hardware.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --sources smooth,sp2_gapped
```

The Slurm job has one Julia process and one CPU thread. Its `process_time.txt`
reports whole-process peak RSS; `samples.csv` contains per-kernel time and
allocation measurements. Inspect `errors.csv` before interpreting any later
rows: failed cases are recorded and the remaining grid continues.
