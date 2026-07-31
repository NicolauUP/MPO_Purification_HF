#!/usr/bin/env bash
# CDW seed-stability test at 256x256:
# three physical points x (+0.1, -0.1, +1e-3) checkerboard seeds.
#SBATCH --job-name=aa-l8-seeds
#SBATCH --output=aa-l8-seeds-%A_%a.out
#SBATCH --error=aa-l8-seeds-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00
#SBATCH --array=1-9%6

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
seed_index=$((task / 3 + 1))

case_labels=(v2_0_u_0p5 v2_0p5_u_0p5 v2_2_u_0p3)
seed_labels=(checkerboard_plus checkerboard_minus checkerboard_weak)
production_moments=(1200 1600 2400)
audit_moments=(1600 2000 3200)
case_label="${case_labels[$((case_index - 1))]}"
seed_label="${seed_labels[$((seed_index - 1))]}"
moments="${production_moments[$((case_index - 1))]}"
final_moments="${audit_moments[$((case_index - 1))]}"

campaign_file="runs/2d/separable_aubry_andre_lside8/campaign_seed_stability.jl"
output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside8_seed_stability/job_${SLURM_ARRAY_JOB_ID}/${case_label}/${seed_label}/M${moments}"
cd "$REPOSITORY_ROOT"
mkdir -p "$(dirname "$output")"

echo "case=$case_label seed=$seed_label M=$moments audit_M=$final_moments"
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

KPM_SEED_STABILITY_CASE="$case_index" \
KPM_SEED_STABILITY_SEED="$seed_index" \
KPM_SCF_MIXER=pulay \
KPM_SCF_PULAY_HISTORY=3 \
KPM_SCF_PULAY_WARMUP=3 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
KPM_SCF_PULAY_COEFFICIENT_LIMIT=16 \
KPM_SCF_PULAY_DAMPING=1.0 \
KPM_SCF_PULAY_DIAGNOSTICS=1 \
KPM_SCF_MOMENTS="$moments" \
KPM_SCF_PROBES=1024 \
KPM_SCF_FINAL_MOMENTS="$final_moments" \
KPM_SCF_FINAL_PROBES=2048 \
KPM_SCF_MAX_ITERATIONS=60 \
  /usr/bin/time -v -o "${output}_process_time.txt" \
  julia --startup-file=no --project=. \
  runs/common/run_square_kpm_scf.jl \
  "$campaign_file" "$output" cuda
