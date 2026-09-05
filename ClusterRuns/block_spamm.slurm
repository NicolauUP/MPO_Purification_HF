#!/bin/bash -l
#SBATCH --job-name=block_spamm
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --mem=32G

set -euo pipefail

: "${N:=1024}"
: "${BLOCK_SIZE:=16}"
: "${V2:=0.5}"
: "${SP2_STEPS:=20}"
: "${TASK:=scf}" # options: diagnostic, sp2, scf
: "${OUTPUT_DIR:=results/block_spamm_${SLURM_JOB_ID:-manual}}"

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPOSITORY_ROOT"
mkdir -p logs "$OUTPUT_DIR"
RUN_OUTPUT="$OUTPUT_DIR/$TASK"

module load julia/1.12.6
module load nvidia-hpc-sdk/24.11
source "$REPOSITORY_ROOT/runs/common/slurm_cuda_runtime.sh"

case "$TASK" in
  diagnostic)
    julia "${JULIA_CUDA_ARGS[@]}" \
      benchmark/sparse_sp2_gpu/diagnose_block_spamm_gpu_frontend.jl \
      "$RUN_OUTPUT" "$N" "$BLOCK_SIZE" "$V2" "$SP2_STEPS"
    ;;
  sp2)
    export SPAMM_PRODUCTION_ONLY=true
    export SPAMM_FULL_SP2_STEPS="$SP2_STEPS"
    julia "${JULIA_CUDA_ARGS[@]}" \
      benchmark/sparse_sp2_gpu/diagnose_block_spamm_gpu_frontend.jl \
      "$RUN_OUTPUT" "$N" "$BLOCK_SIZE" "$V2" 0
    ;;
  scf)
    # This driver fixes BLOCK_SIZE=16 and SP2_STEPS=24 internally.
    julia "${JULIA_CUDA_ARGS[@]}" \
      benchmark/sparse_sp2_gpu/run_block_spamm_scf_1d.jl \
      "$RUN_OUTPUT" "$N" "$V2"
    ;;
  *)
    echo "invalid TASK=$TASK; expected diagnostic, sp2, or scf" >&2
    exit 2
    ;;
esac
