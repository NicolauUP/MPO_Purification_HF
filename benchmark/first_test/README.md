# First cluster test

Submit from the `MPO_Purification_HF` directory:

```bash
sbatch benchmark/first_test/run_tests.slurm
```

The script loads `julia/1.12.6` and runs the maintained
`test/runtests.jl` entry point directly. Direct invocation avoids `Pkg.test()`
trying to download the General registry on the offline compute node. Slurm
expects `JULIA_DEPOT_PATH` to point to the transferred offline depot and
disables automatic environment-wide precompilation. It writes the job output
to:

```text
qtt_hf_first_test-<JOB_ID>.out
qtt_hf_first_test-<JOB_ID>.err
```

This job does not run any benchmark calculation.
