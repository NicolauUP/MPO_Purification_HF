#!/usr/bin/env bash
# Shared runtime setup for MN5 CUDA Julia jobs.
#
# This file is sourced by launchers after loading Julia and the NVIDIA HPC SDK.
# It deliberately makes a pre-built, project-specific Julia depot the only
# depot visible to science jobs. The immutable package/artifact store is linked
# into that depot by submit_cuda_environment_init.slurm.

set -euo pipefail

: "${REPOSITORY_ROOT:?set REPOSITORY_ROOT before sourcing slurm_cuda_runtime.sh}"

export MPO_SHARED_JULIA_DEPOT="${MPO_SHARED_JULIA_DEPOT:-/gpfs/projects/epor78/cluster_depot}"
export MPO_CUDA_RUNTIME_DEPOT="${MPO_CUDA_RUNTIME_DEPOT:-/gpfs/projects/epor78/mpo_julia_runtime_julia1.12_cuda12.6}"

if [[ ! -f "$MPO_CUDA_RUNTIME_DEPOT/.mpohf_cuda_runtime_ready" ]]; then
  echo "CUDA Julia runtime depot is not initialized: $MPO_CUDA_RUNTIME_DEPOT" >&2
  echo "Submit runs/template/submit_cuda_environment_init.slurm once, then retry." >&2
  exit 1
fi

# Science jobs must never write package images or compiled modules. This makes
# an accidental cache miss explicit instead of racing with another job.
export JULIA_DEPOT_PATH="$MPO_CUDA_RUNTIME_DEPOT"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# `nvcc` is installed in the SDK compiler tree, which is not the CUDA toolkit
# root and does not contain NVVM. Resolve the matching CUDA tree from the
# module-provided NVHPC_ROOT instead, then expose exactly that single root.
: "${NVHPC_ROOT:?nvidia-hpc-sdk/24.11 did not define NVHPC_ROOT}"
CUDA_TOOLKIT_VERSION="$(nvcc --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' | tail -n 1)"
test -n "$CUDA_TOOLKIT_VERSION" || { echo "Could not determine CUDA toolkit version." >&2; exit 1; }
CUDA_TOOLKIT_ROOT="$NVHPC_ROOT/cuda/$CUDA_TOOLKIT_VERSION"
test -f "$CUDA_TOOLKIT_ROOT/nvvm/lib64/libnvvm.so" || {
  echo "CUDA toolkit lacks libnvvm: $CUDA_TOOLKIT_ROOT" >&2
  exit 1
}
unset CUDA_PATH CUDA_HOME CUDA_ROOT NVHPC_ROOT
export CUDA_PATH="$CUDA_TOOLKIT_ROOT"
export CUDA_HOME="$CUDA_TOOLKIT_ROOT"
export CUDA_ROOT="$CUDA_TOOLKIT_ROOT"

JULIA_CUDA_ARGS=(--startup-file=no --compiled-modules=existing --pkgimages=existing --project=.)
