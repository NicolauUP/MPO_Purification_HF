#!/usr/bin/env bash
# Shared execution body for the independent H100 SCF scaling submissions.

set -euo pipefail

side_level="${1:?usage: run_scf_scaling_cuda.sh SIDE_LEVEL}"
case "$side_level" in
  7|8|9|10) ;;
  *)
    echo "SIDE_LEVEL must be one of 7, 8, 9, or 10." >&2
    exit 2
    ;;
esac

module load julia/1.12.6
module load nvidia-hpc-sdk/24.11

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0

BENCHMARK_ROOT="${MPO_BENCHMARK_ROOT:?set MPO_BENCHMARK_ROOT to an external result directory}"
OUTPUT_DIR="${BENCHMARK_ROOT}/scf_scaling_lside${side_level}_cuda_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_DIR"

echo "host=$(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "L_side=$side_level L_total=$((2 * side_level))"
echo "SCF scaling output: $OUTPUT_DIR"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

cuda_runtime_version="$(
  nvcc --version |
    sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' |
    tail -n 1
)"
test -n "$cuda_runtime_version" || {
  echo "Could not determine the local CUDA toolkit version from nvcc." >&2
  exit 1
}

MPO_CUDA_RUNTIME_VERSION="$cuda_runtime_version" \
  julia --startup-file=no --project=. -e \
  'using CUDA; CUDA.set_runtime_version!(
      VersionNumber(ENV["MPO_CUDA_RUNTIME_VERSION"]);
      local_toolkit=true,
  )'

/usr/bin/time -v -o "$OUTPUT_DIR/process_time.txt" \
  julia --startup-file=no --project=. \
  benchmark/scf_pilot_2d/scf_pilot_2d.jl \
  --output "$OUTPUT_DIR" \
  --backend cuda \
  --square-fock-method binary_carry \
  --side-level "$side_level" \
  --tx -0.6 \
  --ty -0.35 \
  --U 0.3 \
  --density 0.5 \
  --potential checkerboard \
  --potential-amplitude 0.6 \
  --seed checkerboard_plus \
  --seed-amplitude 0.05 \
  --tci-tol 1e-10 \
  --itensors-tol 1e-10 \
  --maxdim 256 \
  --steps 50 \
  --padding 0.5 \
  --sp2-idempotency-tolerance 2e-4 \
  --sp2-relative-trace-tolerance 1e-6 \
  --observables compact \
  --scf-mixing 0.5 \
  --scf-tol 0.1 \
  --scf-max-iterations 30
