#!/usr/bin/env bash
# Blocked KPM-HF SCF at 1024x512 with safeguarded Pulay mixing.
# Compare its history and audit against the completed linear-mixing 524k run.
#SBATCH --job-name=rect-524k-kpm-pulay
#SBATCH --output=rect-524k-kpm-pulay-%j.out
#SBATCH --error=rect-524k-kpm-pulay-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00

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

output="${MPO_RESULTS_ROOT}/rectangular_1024x512_kpm_scf_pulay/job_${SLURM_JOB_ID}"
cd "$REPOSITORY_ROOT"
echo "host=$(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
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
MPO_CUDA_RUNTIME_VERSION="$cuda_runtime_version" \
  julia --startup-file=no --project=. -e \
  'using CUDA; CUDA.set_runtime_version!(
      VersionNumber(ENV["MPO_CUDA_RUNTIME_VERSION"]);
      local_toolkit=true,
  )'

mkdir -p "$(dirname "$output")"
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

KPM_SCF_MIXER=pulay \
KPM_SCF_PULAY_HISTORY=6 \
KPM_SCF_PULAY_WARMUP=4 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
  /usr/bin/time -v -o "${output}_process_time.txt" \
  julia --startup-file=no --project=. \
  runs/common/run_rectangular_kpm_scf.jl "$output" cuda
