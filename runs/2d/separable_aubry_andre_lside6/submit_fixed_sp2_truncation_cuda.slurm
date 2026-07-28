#!/usr/bin/env bash
# Tight-bound, chi=512 SP2 trajectory diagnostic. Array tasks compare the
# operational truncation cutoffs 1e-8 and 1e-10 on independent H100s.
#SBATCH --job-name=qtt_aa_l6_sp2_trunc
#SBATCH --output=aa-l6-sp2-trunc-%A_%a.out
#SBATCH --error=aa-l6-sp2-trunc-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-2

set -euo pipefail
module load julia/1.12.6
module load nvidia-hpc-sdk/24.11

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
NVCC_PATH="$(readlink -f "$(command -v nvcc)")"
CUDA_TOOLKIT_ROOT="$(dirname "$(dirname "$NVCC_PATH")")"
export CUDA_PATH="$CUDA_TOOLKIT_ROOT"
export CUDA_HOME="$CUDA_TOOLKIT_ROOT"
export CUDA_ROOT="$CUDA_TOOLKIT_ROOT"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=16
export MKL_NUM_THREADS=16
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

CUTOFFS=(1e-8 1e-10)
LABELS=(1e-8 1e-10)
cutoff="${CUTOFFS[$((SLURM_ARRAY_TASK_ID - 1))]}"
label="${LABELS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside6_sp2_truncation/cutoff_${label}"

cd "$REPOSITORY_ROOT"
mkdir -p "$(dirname "$output")"
echo "host=$(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "CUDA_PATH=$CUDA_PATH"
echo "fixed-H SP2 maxdim=512 cutoff=$cutoff bounds=[-4.1,4.1]"
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

memory_log="$(dirname "$output")/gpu_memory_cutoff_${label}.csv"
echo '"timestamp","memory_used_mib","memory_free_mib"' > "$memory_log"
nvidia-smi \
  --query-gpu=timestamp,memory.used,memory.free \
  --format=csv,noheader,nounits \
  --loop=1 >> "$memory_log" &
sampler_pid=$!
cleanup_sampler() {
  kill "$sampler_pid" 2>/dev/null || true
  wait "$sampler_pid" 2>/dev/null || true
}
trap cleanup_sampler EXIT

/usr/bin/time -v -o "$(dirname "$output")/process_time_cutoff_${label}.txt" \
  julia --startup-file=no --project=. \
  runs/common/compare_square_sp2_dense_projector.jl \
  runs/2d/separable_aubry_andre_lside6/campaign_tight_fixed_h.jl \
  1 512 "$output" 2e-4 cuda "$cutoff" 1e-6
