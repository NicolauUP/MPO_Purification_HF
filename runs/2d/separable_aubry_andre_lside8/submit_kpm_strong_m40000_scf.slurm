#!/usr/bin/env bash
# Full high-degree KPM Hartree--Fock SCF for the strong quasiperiodic CDW.
# Production: M=40000, R=2048. Independent final audit: M=60000, R=2048.
# The fixed-field degree ladder established that R=2048 is not the limiting
# error source, so no R=4096 production pass is requested here.
#SBATCH --job-name=aa-l8-m40k-scf
#SBATCH --output=aa-l8-m40k-scf-%j.out
#SBATCH --error=aa-l8-m40k-scf-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=08:00:00

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

output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside8_strong_m40000_scf/job_${SLURM_JOB_ID}"
cd "$REPOSITORY_ROOT"
mkdir -p "$(dirname "$output")"
echo "output=$output"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader

NVCC_PATH="$(readlink -f "$(command -v nvcc)")"
CUDA_TOOLKIT_ROOT="$(dirname "$(dirname "$NVCC_PATH")")"
export CUDA_PATH="$CUDA_TOOLKIT_ROOT"
export CUDA_HOME="$CUDA_TOOLKIT_ROOT"
export CUDA_ROOT="$CUDA_TOOLKIT_ROOT"
cuda_runtime_version="$(nvcc --version | sed -n 's/.*release \([0-9][0-9]*\.[0-9][0-9]*\).*/\1/p' | tail -n 1)"
test -n "$cuda_runtime_version" || { echo "Could not determine CUDA runtime" >&2; exit 1; }
MPO_CUDA_RUNTIME_VERSION="$cuda_runtime_version" julia --startup-file=no --project=. -e \
  'using CUDA; CUDA.set_runtime_version!(VersionNumber(ENV["MPO_CUDA_RUNTIME_VERSION"]); local_toolkit=true)'

memory_log="${output}_gpu_memory.csv"
echo '"timestamp","memory_used_mib","memory_free_mib"' > "$memory_log"
nvidia-smi --query-gpu=timestamp,memory.used,memory.free --format=csv,noheader,nounits --loop=1 >> "$memory_log" &
sampler_pid=$!
trap 'kill "$sampler_pid" 2>/dev/null || true; wait "$sampler_pid" 2>/dev/null || true' EXIT

KPM_SCF_MIXER=pulay \
KPM_SCF_PULAY_HISTORY=3 \
KPM_SCF_PULAY_WARMUP=3 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
KPM_SCF_PULAY_COEFFICIENT_LIMIT=16 \
KPM_SCF_PULAY_DAMPING=1.0 \
KPM_SCF_PULAY_DIAGNOSTICS=1 \
KPM_SCF_MOMENTS=40000 \
KPM_SCF_PROBES=2048 \
KPM_SCF_FINAL_MOMENTS=60000 \
KPM_SCF_FINAL_PROBES=2048 \
KPM_SCF_MAX_ITERATIONS=60 \
  /usr/bin/time -v -o "${output}_process_time.txt" \
  julia --startup-file=no --project=. \
  runs/common/run_square_kpm_scf.jl \
  runs/2d/separable_aubry_andre_lside8/campaign_strong_m40000.jl \
  "$output" cuda
