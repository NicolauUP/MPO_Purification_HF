#!/bin/bash -l
#SBATCH --job-name=cusparse_csr
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=20
#SBATCH --time=01:00:00
#SBATCH --account=epor78
#SBATCH --qos=acc_ehpc

set -euo pipefail

: "${N:=1024}"
: "${V2:=2.0}"
: "${U:=0.3}"
: "${SEED:=0.5}"
: "${CUTOFF:=1e-8}"
: "${TASK:=scf}" # options: scale, scf
: "${OUTPUT_DIR:=results/cusparse_${SLURM_JOB_ID:-manual}}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/../Project.toml" ]]; then
  export REPOSITORY_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/Project.toml" ]]; then
  export REPOSITORY_ROOT="${SLURM_SUBMIT_DIR}"
else
  export REPOSITORY_ROOT="${REPOSITORY_ROOT:-$(pwd)}"
fi
cd "${REPOSITORY_ROOT}"
mkdir -p logs "$OUTPUT_DIR"
RUN_OUTPUT="$OUTPUT_DIR/$TASK"

module load julia/1.12.6
module load nvidia-hpc-sdk/24.11
source "$REPOSITORY_ROOT/runs/common/slurm_cuda_runtime.sh"

case "$TASK" in
  scale)
    julia "${JULIA_CUDA_ARGS[@]}" \
      benchmark/sparse_sp2_gpu/scale_sparse_sp2_1d.jl \
      "$RUN_OUTPUT" "$N" "$V2" "$CUTOFF"
    ;;
  scf)
    julia "${JULIA_CUDA_ARGS[@]}" \
      benchmark/sparse_sp2_gpu/run_sparse_sp2_scf_1d.jl \
      "$RUN_OUTPUT" "$N" "$V2" "$U" "$SEED" "$CUTOFF"
    ;;
  *)
    echo "invalid TASK=$TASK; expected scale or scf" >&2
    exit 2
    ;;
esac
