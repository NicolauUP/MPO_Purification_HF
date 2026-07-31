#!/usr/bin/env bash
# Fixed-resolution KPM-SCF ladder at 256x256 for three representative points.
# 3 physical points x 4 moment values x 2 initial seeds = 24 tasks.
#SBATCH --job-name=aa-l8-moments
#SBATCH --output=aa-l8-moments-%A_%a.out
#SBATCH --error=aa-l8-moments-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00
#SBATCH --array=1-24%6

set -euo pipefail
module load julia/1.12.6
module load nvidia-hpc-sdk/24.11

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root}"

task=$((SLURM_ARRAY_TASK_ID - 1))
case_index=$((task % 3 + 1))
moment_index=$((task / 3 % 4))
seed_index=$((task / 12))
moments=(400 800 1200 1600)
production_moments="${moments[$moment_index]}"

case "$case_index" in
  1) case_label="v2_0_u_0p5" ;;
  2) case_label="v2_0p5_u_0p5" ;;
  3) case_label="v2_2_u_0p3" ;;
esac
case "$seed_index" in
  0) seed_mode="checkerboard" ;;
  1) seed_mode="uniform" ;;
  *) echo "invalid seed index: $seed_index" >&2; exit 2 ;;
esac

campaign_file="runs/2d/separable_aubry_andre_lside8/campaign_moment_convergence.jl"
output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside8_moment_convergence/job_${SLURM_ARRAY_JOB_ID}/${case_label}/${seed_mode}/M${production_moments}"
cd "$REPOSITORY_ROOT"
mkdir -p "$(dirname "$output")"

echo "case=$case_label case_index=$case_index seed=$seed_mode moments=$production_moments"
echo "output=$output"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

NVCC_PATH="$(readlink -f "$(command -v nvcc)")"
CUDA_TOOLKIT_ROOT="$(dirname "$(dirname "$NVCC_PATH")")"
export CUDA_PATH="$CUDA_TOOLKIT_ROOT"
export CUDA_HOME="$CUDA_TOOLKIT_ROOT"
export CUDA_ROOT="$CUDA_TOOLKIT_ROOT"
cuda_runtime_version="$(
  nvcc --version |
    sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' |
    tail -n 1
)"
test -n "$cuda_runtime_version" || {
  echo "Could not determine local CUDA toolkit version." >&2
  exit 1
}
MPO_CUDA_RUNTIME_VERSION="$cuda_runtime_version" \
  julia --startup-file=no --project=. -e \
  'using CUDA; CUDA.set_runtime_version!(
      VersionNumber(ENV["MPO_CUDA_RUNTIME_VERSION"]);
      local_toolkit=true,
  )'

memory_log="${output}_gpu_memory.csv"
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

KPM_MOMENT_CASE="$case_index" \
KPM_MOMENT_SEED="$seed_mode" \
KPM_SCF_MIXER=pulay \
KPM_SCF_PULAY_HISTORY=3 \
KPM_SCF_PULAY_WARMUP=3 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
KPM_SCF_PULAY_COEFFICIENT_LIMIT=16 \
KPM_SCF_PULAY_DAMPING=1.0 \
KPM_SCF_MOMENTS="$production_moments" \
KPM_SCF_PROBES=1024 \
KPM_SCF_FINAL_MOMENTS=2000 \
KPM_SCF_FINAL_PROBES=2048 \
  /usr/bin/time -v -o "${output}_process_time.txt" \
  julia --startup-file=no --project=. \
  runs/common/run_square_kpm_scf.jl \
  "$campaign_file" "$output" cuda
