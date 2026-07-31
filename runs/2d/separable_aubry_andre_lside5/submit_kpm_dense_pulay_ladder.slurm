#!/usr/bin/env bash
# Three safeguarded Pulay/DIIS settings, each cross-validated with dense HF.
#SBATCH --job-name=aa-l5-kpm-pulay
#SBATCH --output=aa-l5-kpm-pulay-%A_%a.out
#SBATCH --error=aa-l5-kpm-pulay-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=02:00:00
#SBATCH --array=1-3

set -euo pipefail
module load julia/1.12.6
module load nvidia-hpc-sdk/24.11

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-20}"
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root}"

case "${SLURM_ARRAY_TASK_ID}" in
  1) label="history4_cap16"; history=4; coefficient_limit=16 ;;
  2) label="history4_cap32"; history=4; coefficient_limit=32 ;;
  3) label="history3_cap16"; history=3; coefficient_limit=16 ;;
  *) echo "unknown array task: ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

root="${MPO_RESULTS_ROOT}/kpm_dense_lside5_pulay_ladder/job_${SLURM_ARRAY_JOB_ID}/${label}"
kpm_output="$root/kpm"
campaign_file="runs/2d/separable_aubry_andre_lside5/campaign_v2_0p5_strict.jl"
cd "$REPOSITORY_ROOT"
mkdir -p "$root"
echo "label=$label history=$history coefficient_limit=$coefficient_limit"

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

KPM_SCF_MIXER=pulay \
KPM_SCF_PULAY_DAMPING=1.0 \
KPM_SCF_PULAY_HISTORY="$history" \
KPM_SCF_PULAY_WARMUP=4 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
KPM_SCF_PULAY_COEFFICIENT_LIMIT="$coefficient_limit" \
KPM_SCF_PULAY_DIAGNOSTICS=1 \
KPM_SCF_MOMENTS=1200 \
KPM_SCF_PROBES=1024 \
KPM_SCF_FINAL_MOMENTS=1600 \
KPM_SCF_FINAL_PROBES=1024 \
  julia --startup-file=no --project=. \
  runs/common/run_square_kpm_scf.jl \
  "$campaign_file" "$kpm_output" cuda

dense_root="$root/dense_root"
MPO_RESULTS_ROOT="$dense_root" \
  julia --startup-file=no --project=. \
  runs/common/execute_square_dense_hf_campaign.jl "$campaign_file" 1

dense_output="$dense_root/separable_aubry_andre_lside5_v2_0p5_strict_dense_hf/task_0001_v2_0p5_u_1"
julia --startup-file=no --project=. \
  runs/common/compare_kpm_dense_outputs.jl \
  "$kpm_output" "$dense_output" "$root/comparison.toml"

cat "$root/comparison.toml"
