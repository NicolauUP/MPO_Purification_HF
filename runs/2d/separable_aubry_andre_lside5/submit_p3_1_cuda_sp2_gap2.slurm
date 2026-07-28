#!/usr/bin/env bash
# P3.1: first CUDA correctness/timing calibration. This repeats the validated
# fixed-Hamiltonian gap-2 SP2 case on one MN5 ACC GPU.
#SBATCH --job-name=qtt_sp2_cuda_gap2
#SBATCH --output=p3-1-sp2-cuda-gap2-%j.out
#SBATCH --error=p3-1-sp2-cuda-gap2-%j.err
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
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

output="${MPO_RESULTS_ROOT}/square_lside5_sp2_cuda_gap2/maxdim_256_tol_1e-10"

cd "$REPOSITORY_ROOT"
echo "host=$(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
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
echo "Configuring CUDA.jl to use local CUDA ${cuda_runtime_version}"

MPO_CUDA_RUNTIME_VERSION="$cuda_runtime_version" \
  julia --startup-file=no --project=. -e \
  'using CUDA; CUDA.set_runtime_version!(
      VersionNumber(ENV["MPO_CUDA_RUNTIME_VERSION"]);
      local_toolkit=true,
  )'

/usr/bin/time -v julia --startup-file=no --project=. \
  runs/common/compare_square_sp2_dense_projector.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  2 256 "$output" 2e-4 cuda 1e-10 1e-6
