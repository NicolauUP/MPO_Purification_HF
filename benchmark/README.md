# Extraction-scaling benchmark

`extraction_scaling.jl` compares the default cached-TCI and experimental
binary-carry 1D extractors for both Hartree and Fock fields. It is a
performance/equivalence benchmark, not a self-consistent production
calculation.

## Cluster preflight

First run a small case to confirm the cluster environment and output path:

```bash
cd /gpfs/projects/epor78/MPO_HartreeFock_Purification/MPO_Purification_HF
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
mkdir -p "$MPO_BENCHMARK_ROOT"

bash benchmark/run_extraction_scaling.sh \
  "$MPO_BENCHMARK_ROOT/preflight" \
  --levels 4,6 --sources smooth --warmups 1 \
  --repetitions-small 1 --repetitions-large 1
```

This is intentionally a direct command, not a full Slurm benchmark. It
creates a small result directory and verifies that all four extractors can run.

## Full cluster benchmark

The Slurm wrapper requires an external result root and submits the default
grid: `L=4,6,8,10,12,14`, sources `smooth,sp2_gapped`, one warmup per
extractor, five repetitions for `L≤10`, and three for `L≥12`.

```bash
export MPO_BENCHMARK_ROOT=/gpfs/projects/epor78/MPO_HF_benchmarks
sbatch benchmark/extraction_scaling.slurm
```

This performs 208 timed extractor calls plus 48 warmups. The `L=14` cases are
the limiting memory/time points. The wrapper deliberately uses one Julia
thread because it measures sequential extractor implementations; do not
interpret it as a many-core scaling study.

## Output

Each full job writes `$MPO_BENCHMARK_ROOT/extraction_<job-id>/` containing:

- `metadata.txt` — host, Julia, Git, benchmark configuration;
- `samples.csv` — raw timing/allocation samples;
- `summary.csv` — per-size/source/extractor medians and equivalence metrics;
- `errors.csv` — case-level failures, while later cases continue;
- `README.md` — output-field description;
- `process_time.txt` — whole-process OS memory/time report.

`process_time.txt` is a process-wide peak RSS, not a per-extractor memory
measurement. Use `summary.csv` and `samples.csv` for method-level comparison.
