# CUDA Slurm submission policy

Use `submit_mpo_cuda.slurm` as the only starting point for new MPO CUDA
campaign submissions on MareNostrum5.

The template was distilled from the cluster's successful public MPO and C8
campaign launchers. It intentionally fixes the runtime details that have
otherwise drifted between scripts:

- load Julia 1.12.6 and NVIDIA HPC SDK 24.11;
- use an isolated, pre-built runtime depot whose packages/artifacts are
  symlinked from the shared offline depot;
- derive the actual CUDA toolkit (including NVVM) from module-provided
  `NVHPC_ROOT` and the `nvcc` version, then replace every CUDA root variable
  with that single toolkit root;
- require the one-off `submit_cuda_environment_init.slurm` job to configure
  CUDA.jl and build package caches before any science job;
- run science with existing compiled modules and pkgimages only, so it cannot
  mutate or race over the shared cache.

Do not copy old submitters, run `CUDA.set_runtime_version!` in a science job,
introduce another Julia depot, or use `--compiled-modules=no`. If the runtime
cache is missing or must be rebuilt after a dependency/runtime change, submit
the initializer alone and wait for it to succeed. Place physics and numerical
parameters in a campaign source; configure a new run through the `MPOHF_*`
environment variables at submission time.

For a nonstandard resource request, make a thin wrapper that changes only the
`#SBATCH` allocation lines and invokes this template. It must not duplicate
the CUDA or Julia setup.
