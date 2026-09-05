#!/usr/bin/env bash
#SBATCH --job-name=square_diag
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=20
#SBATCH --time=01:00:00
#SBATCH --account=epor78
#SBATCH --qos=acc_ehpc
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

: "${TASK:=fixed_sp2}"       # options: fixed_sp2, mcweeny, sp2_dense, projector_compression
: "${OUTPUT_DIR:=results/diag_${SLURM_JOB_ID:-manual}}"
: "${CAMPAIGN:=runs/template/campaign.jl}"
: "${CAMPAIGN_TASK:=1}"
: "${MAXDIM:=256}"
: "${BACKEND:=cuda}"
: "${SPECTRAL_LOWER:=-4.25}"
: "${SPECTRAL_UPPER:=4.25}"
: "${MU:=0.0}"
: "${ITENSORS_TOL:=1e-10}"
: "${RECURRENCE:=standard}"
: "${SP2_IDEMPOTENCY_TOLERANCE:=1e-3}"
: "${SP2_RELATIVE_TRACE_TOLERANCE:=1e-6}"
: "${PROJECTOR_MODE:=exact}"
: "${CUTOFFS:=1e-14,1e-12,1e-10,1e-8,1e-6}"

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/../Project.toml" ]]; then
  export REPOSITORY_ROOT="$(cd "${SLURM_SUBMIT_DIR}/.." && pwd)"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]] && [[ -f "${SLURM_SUBMIT_DIR}/Project.toml" ]]; then
  export REPOSITORY_ROOT="${SLURM_SUBMIT_DIR}"
else
  export REPOSITORY_ROOT="${REPOSITORY_ROOT:-$(pwd)}"
fi
cd "${REPOSITORY_ROOT}"
mkdir -p logs "$OUTPUT_DIR"
[ -f runs/common/slurm_cuda_runtime.sh ] && source runs/common/slurm_cuda_runtime.sh

case "$TASK" in
  fixed_sp2)
    julia --project=. runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
      "$CAMPAIGN" "$CAMPAIGN_TASK" "$MAXDIM" "$OUTPUT_DIR/fixed_sp2" \
      "$BACKEND" "$SPECTRAL_LOWER" "$SPECTRAL_UPPER"
    ;;
  mcweeny)
    julia --project=. runs/common/diagnose_square_mcweeny_mu_fixed_hamiltonian.jl \
      "$CAMPAIGN" "$CAMPAIGN_TASK" "$MAXDIM" "$OUTPUT_DIR/mcweeny" \
      "$SPECTRAL_LOWER" "$SPECTRAL_UPPER" "$MU" "$ITENSORS_TOL" "$RECURRENCE"
    ;;
  sp2_dense)
    julia --project=. runs/common/compare_square_sp2_dense_projector.jl \
      "$CAMPAIGN" "$CAMPAIGN_TASK" "$MAXDIM" "$OUTPUT_DIR/sp2_dense" \
      "$SP2_IDEMPOTENCY_TOLERANCE" "$BACKEND" "$ITENSORS_TOL" \
      "$SP2_RELATIVE_TRACE_TOLERANCE"
    ;;
  projector_compression)
    julia --project=. runs/common/compare_square_projector_compression.jl \
      "$CAMPAIGN" "$CAMPAIGN_TASK" "$PROJECTOR_MODE" "$MAXDIM" \
      "$OUTPUT_DIR/projector_compression" "$CUTOFFS" "$ITENSORS_TOL" \
      "$SP2_IDEMPOTENCY_TOLERANCE" "$SP2_RELATIVE_TRACE_TOLERANCE"
    ;;
  *)
    echo "error: TASK must be fixed_sp2, mcweeny, sp2_dense, or projector_compression (got: $TASK)" >&2
    exit 2
    ;;
esac
