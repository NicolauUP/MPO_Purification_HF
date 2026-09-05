#!/usr/bin/env bash
#SBATCH --job-name=mpohf_cuda
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

: "${CAMPAIGN:=runs/template/campaign.jl}"
: "${METHOD:=mpo}"           # options: mpo, dense, kpm
: "${TASK:=1}"
: "${OUTPUT_ROOT:=results/mpohf_${SLURM_JOB_ID:-manual}}"
: "${DEFER_OBSERVABLES:=false}"
: "${VERBOSE:=true}"

mkdir -p logs "$OUTPUT_ROOT"
[ -f runs/common/slurm_cuda_runtime.sh ] && source runs/common/slurm_cuda_runtime.sh

julia --project=. bin/mpohf.jl --campaign "$CAMPAIGN" --task "$TASK" --method "$METHOD" --output "$OUTPUT_ROOT"
