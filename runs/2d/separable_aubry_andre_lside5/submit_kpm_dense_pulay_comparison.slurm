#!/usr/bin/env bash
# Safeguarded Pulay KPM-SCF versus independent dense-HF comparison at N=1024.
# This is the validation gate before using Pulay on large blocked KPM runs.
#SBATCH --job-name=aa-l5-kpm-pulay
#SBATCH --output=aa-l5-kpm-pulay-%j.out
#SBATCH --error=aa-l5-kpm-pulay-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=02:00:00

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

root="${MPO_RESULTS_ROOT}/kpm_dense_lside5_pulay/job_${SLURM_JOB_ID}"
kpm_output="$root/kpm"
campaign_file="runs/2d/separable_aubry_andre_lside5/campaign_v2_0p5_strict.jl"
cd "$REPOSITORY_ROOT"
mkdir -p "$root"

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
KPM_SCF_PULAY_HISTORY=6 \
KPM_SCF_PULAY_WARMUP=4 \
KPM_SCF_PULAY_REGULARIZATION=1e-12 \
KPM_SCF_MOMENTS=1200 \
KPM_SCF_PROBES=1024 \
KPM_SCF_FINAL_MOMENTS=1600 \
KPM_SCF_FINAL_PROBES=1024 \
  julia --startup-file=no --project=. \
  runs/common/run_square_kpm_scf.jl \
  "$campaign_file" "$kpm_output" cuda

comparison_results_root="$root/dense_root"
MPO_RESULTS_ROOT="$comparison_results_root" \
  julia --startup-file=no --project=. \
  runs/common/execute_square_dense_hf_campaign.jl "$campaign_file" 1

dense_output="$comparison_results_root/separable_aubry_andre_lside5_v2_0p5_strict_dense_hf/task_0001_v2_0p5_u_1"
julia --startup-file=no --project=. \
  runs/common/compare_kpm_dense_outputs.jl \
  "$kpm_output" "$dense_output" "$root/comparison.toml"

cat "$root/comparison.toml"
