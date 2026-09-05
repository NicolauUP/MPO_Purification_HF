#!/usr/bin/env bash
#SBATCH --job-name=kpm_scf
#SBATCH --gres=gpu:1
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

: "${GEOMETRY:=square}"      # options: square, rectangular
: "${OUTPUT_DIR:=results/kpm_${SLURM_JOB_ID:-manual}}"
: "${CAMPAIGN:=runs/template/campaign.jl}"
: "${BACKEND:=cuda}"

mkdir -p logs "$OUTPUT_DIR"
[ -f runs/common/slurm_cuda_runtime.sh ] && source runs/common/slurm_cuda_runtime.sh

case "$GEOMETRY" in
  square)
    julia --project=. runs/common/run_square_kpm_scf.jl \
      "$CAMPAIGN" "$OUTPUT_DIR/square" "$BACKEND"
    ;;
  rectangular)
    julia --project=. runs/common/run_rectangular_kpm_scf.jl \
      "$OUTPUT_DIR/rectangular" "$BACKEND"
    ;;
  *)
    echo "error: GEOMETRY must be square or rectangular (got: $GEOMETRY)" >&2
    exit 2
    ;;
esac
