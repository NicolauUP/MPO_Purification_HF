#!/usr/bin/env bash
#SBATCH --job-name=mpohf_cuda
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=20
#SBATCH --time=02:00:00
#SBATCH --account=epor78
#SBATCH --qos=acc_ehpc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

: "${CAMPAIGN:=runs/template/campaign.jl}"
: "${METHOD:=mpo}"           # options: mpo, dense, kpm
: "${TASK:=1}"
: "${OUTPUT_ROOT:=results/mpohf_${SLURM_JOB_ID:-manual}}"
: "${DEFER_OBSERVABLES:=false}"
: "${VERBOSE:=true}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/../Project.toml" ]]; then
  export REPOSITORY_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/Project.toml" ]]; then
  export REPOSITORY_ROOT="${SLURM_SUBMIT_DIR}"
else
  export REPOSITORY_ROOT="${REPOSITORY_ROOT:-$(pwd)}"
fi
cd "${REPOSITORY_ROOT}"
mkdir -p logs "$OUTPUT_ROOT"
[ -f runs/common/slurm_cuda_runtime.sh ] && source runs/common/slurm_cuda_runtime.sh

julia --project=. bin/mpohf.jl --campaign "$CAMPAIGN" --task "$TASK" --method "$METHOD" --output "$OUTPUT_ROOT"
