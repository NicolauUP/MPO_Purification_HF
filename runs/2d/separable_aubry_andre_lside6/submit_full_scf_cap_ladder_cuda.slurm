#!/usr/bin/env bash
# Full 64x64 quasiperiodic HF calculation at three MPO caps. Run this after
# inspecting the fixed-H cap ladder; task 3 may exceed one H100's memory.
#SBATCH --job-name=qtt_aa_l6_scf_caps
#SBATCH --output=aa-l6-scf-caps-%A_%a.out
#SBATCH --error=aa-l6-scf-caps-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-3

set -euo pipefail
module load julia/1.12.6
module load nvidia-hpc-sdk/24.11

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_BENCHMARK_ROOT="${MPO_BENCHMARK_ROOT:?set an external benchmark root before submitting}"

MAXDIMS=(256 512 1024)
maxdim="${MAXDIMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_BENCHMARK_ROOT}/separable_aubry_andre_lside6_scf/maxdim_${maxdim}"

cd "$REPOSITORY_ROOT"
ispath_message="refusing to overwrite existing output: $output"
test ! -e "$output" || { echo "$ispath_message" >&2; exit 1; }
mkdir -p "$output"
echo "host=$(hostname)"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}"
echo "full SCF maxdim=$maxdim output=$output"
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

memory_log="$output/gpu_memory_timeseries.csv"
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

/usr/bin/time -v -o "$output/process_time.txt" \
  julia --startup-file=no --project=. \
  benchmark/scf_pilot_2d/scf_pilot_2d.jl \
  --output "$output" \
  --backend cuda \
  --square-fock-method binary_carry \
  --side-level 6 \
  --hopping separable_aubry_andre \
  --hopping-modulation 0.5 \
  --hopping-frequency 0.5808802290397618 \
  --U 1.0 \
  --density 0.5 \
  --potential none \
  --seed checkerboard_plus \
  --seed-amplitude 0.1 \
  --tci-tol 1e-10 \
  --itensors-tol 1e-10 \
  --maxdim "$maxdim" \
  --steps 50 \
  --padding 0.5 \
  --spectral-lower -8.5 \
  --spectral-upper 12.5 \
  --sp2-idempotency-tolerance 2e-4 \
  --sp2-relative-trace-tolerance 1e-6 \
  --observables full \
  --scf-mixing 0.5 \
  --scf-tol 0.1 \
  --scf-max-iterations 30
