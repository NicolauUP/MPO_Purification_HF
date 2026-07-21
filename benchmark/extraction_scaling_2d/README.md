# Square-lattice extraction benchmark

This benchmark measures two open-square **binary-carry Hartree** extractors on
`2^L_side × 2^L_side` lattices. Thus `L_side=10` means `L_total=20` QTT bits
and `N=1,048,576` sites.

It is deliberately not a full SCF calculation. Density preparation happens
outside the timed region. The timed kernels are:

- `hartree_four_carry` — four direct interleaved-QTT carry shifts followed by
  a compressed MPO sum.
- `hartree_adjacency` — one fused QTT adjacency MPO implementing
  `U * A_nn * n`, with direction and carry stored in its virtual state. It
  avoids that four-MPO sum.

For every result, the benchmark measures the Hartree value directly from the
density MPO at four corners and one bulk site. `summary.csv` records the
maximum and mean absolute error of those probes. This is the accuracy
diagnostic; it does not rely on agreement with another reconstructed MPO.

There is currently no square binary-carry Fock implementation. This benchmark
therefore makes no Fock performance or accuracy claim.

## Cluster preflight

Submit the preflight first. It runs an exact checkerboard density on a `4×4`
lattice and makes no large allocation. This source has
`n(x,y)=0.5+0.2*(-1)^(x+y)` and does not depend on 2D SP2.

```bash
cd /gpfs/projects/epor78/MPO_HartreeFock_Purification/MPO_Purification_HF
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
mkdir -p "$MPO_BENCHMARK_ROOT"

sbatch benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --side-levels 2 --sources checkerboard_exact --warmups 1 \
  --repetitions-small 1 --repetitions-large 1
```

## Full scaling run

The default grid is `L_side=2,3,...,10`, equivalent to total QTT bit counts
`L_total=4,6,...,20`. `checkerboard_exact` is the primary controlled source:
it is bounded, CDW-like, compact in QTT, and independent of purification.
Use it to compare the two extraction methods directly. The `smooth` source
remains useful as a second structured density.

`sp2_gapped` prepares a static checkerboard CDW-potential density,
`W(x,y)=+0.6` on one sublattice and `-0.6` on the other. It is SP2 only, not
an SCF calculation. Current square-MPO SP2 may stagnate, so it is a separate
purification feasibility diagnostic—not the primary Hartree benchmark.

```bash
sbatch --export=ALL,MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks \
  benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --side-levels 4,6,8 --sources sp2_gapped \
  --warmups 0 --repetitions-small 1 --repetitions-large 1
```

This performs SP2 plus both binary-carry extractions at each size, without SCF.
Inspect `errors.csv`, direct-probe errors, and bond dimensions before advancing
to `L_side=10`.

To diagnose the known SP2 cutoff sensitivity at the first cluster-sized
checkerboard case (`L_side=4`, i.e. a `16 x 16` lattice), pass the MPO cutoff
explicitly:

```bash
sbatch --export=ALL,MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks \
  benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --side-levels 4 --sources sp2_gapped --itensors-tol 1e-14 \
  --warmups 0 --repetitions-small 1 --repetitions-large 1
```

This is a diagnosis, not a production SCF run. Compare its `errors.csv` and
SP2 source details with the corresponding `1e-12` run; do not infer a general
production cutoff from one lattice size.

## Exact dense validation of the cluster checkerboard case

`sp2_dense_validation_2d.jl` is a deliberately bounded diagnostic. It allows
only `L_side <= 4`; its largest dense matrix is therefore the `256 x 256`
checkerboard Hamiltonian. It uses dense diagonalization to validate the Fermi
gap and spectral bounds, then compares the final MPO-SP2 density against the
exact occupied projector. It writes `summary.csv`, `sp2_history.txt`,
`metadata.txt`, and `process_time.txt`.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
  benchmark/extraction_scaling_2d/sp2_dense_validation_2d.slurm \
  --side-level 4 --itensors-tol 1e-14 --itensors-maxdim 128 \
  --steps 50 --padding 0.5
```

This job is the appropriate way to decide whether the apparently converged
`16 x 16` MPO-SP2 density is also close to the dense projector. Do not raise
`--side-level` above 4: the script rejects that request to avoid an accidental
large dense diagonalization.

## Refinement after the current SP2 stopping rule

`sp2_refinement_validation_2d.jl` is a separate diagnostic. It first runs the
unchanged production SP2 procedure, then applies extra polynomials only inside
the diagnostic and writes dense-reference metrics after each one. It tests
whether the hard-coded SP2 stopping threshold returns too early.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
  benchmark/extraction_scaling_2d/sp2_refinement_validation_2d.slurm \
  --side-level 4 --itensors-tol 1e-14 --itensors-maxdim 256 \
  --steps 50 --padding 0.5 --extra-iterations 8
```

Read `refinement.csv`: row zero is the normal SP2 stopping state, and rows
1--8 are forced extra polynomials. Production source files and normal SP2
behavior are unchanged by this diagnostic.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch benchmark/extraction_scaling_2d/extraction_scaling_2d.slurm \
  --sources checkerboard_exact,smooth
```

The Slurm job has one Julia process and one CPU thread. Its `process_time.txt`
reports whole-process peak RSS; `samples.csv` contains per-kernel time and
allocation measurements. Inspect `errors.csv` and the direct-probe error
columns before interpreting any performance result: failed cases are recorded
and the remaining grid continues.

## Explicit adjacency-compression sweep

`adjacency_compression_sweep.jl` starts from the accurate, untruncated fused
adjacency field for the smooth synthetic density. It then compresses copies of
that field under explicit `(cutoff, maxdim)` policies. Its `summary.csv` is
plot-ready: plot direct-probe error against `field_max_chi`,
`median_allocations_bytes`, or `median_time_s`.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch --export=ALL,MPO_BENCHMARK_ROOT="$MPO_BENCHMARK_ROOT" \
  benchmark/extraction_scaling_2d/adjacency_compression_sweep.slurm \
  --side-levels 4,6,8,10 \
  --cutoffs 1e-14,1e-12,1e-10,1e-8 \
  --maxdims 16,32,64,128,256 \
  --repetitions 3
```

Record each completed run and its interpretation in the tracked
[`docs/logbook/benchmark_runs.md`](../../../docs/logbook/benchmark_runs.md).
