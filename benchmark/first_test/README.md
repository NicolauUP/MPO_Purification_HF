# First cluster test

Submit from the `MPO_Purification_HF` directory:

```bash
sbatch benchmark/first_test/run_tests.slurm
```

The script loads `julia/1.12.6` and runs the maintained package test suite with
`Pkg.test()`. Slurm writes the job output to:

```text
qtt_hf_first_test-<JOB_ID>.out
qtt_hf_first_test-<JOB_ID>.err
```

This job does not run any benchmark calculation.
