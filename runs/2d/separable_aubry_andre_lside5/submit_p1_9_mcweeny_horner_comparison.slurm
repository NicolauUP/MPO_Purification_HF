#!/usr/bin/env bash
# P1.9: standard versus Horner-form fixed-mu McWeeny evaluation.
# Tasks 1/2 use gap 0.2; tasks 3/4 use gap 1.0. Each pair differs only in
# whether 3P^2-2P^3 is evaluated directly or as P^2(3I-2P).
#SBATCH --job-name=qtt_aa_2d_p19
#SBATCH --output=p1-9-horner-%A_%a.out
#SBATCH --error=p1-9-horner-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-4

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

CAMPAIGN_TASKS=(1 1 2 2)
GAP_LABELS=(gap_0p2 gap_0p2 gap_1p0 gap_1p0)
RECURRENCES=(standard horner standard horner)
index=$((SLURM_ARRAY_TASK_ID - 1))
campaign_task="${CAMPAIGN_TASKS[$index]}"
gap_label="${GAP_LABELS[$index]}"
recurrence="${RECURRENCES[$index]}"
output="${MPO_RESULTS_ROOT}/square_lside5_mcweeny_horner_comparison/${gap_label}/${recurrence}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/diagnose_square_mcweeny_mu_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_horner_comparison.jl \
  "$campaign_task" 256 "$output" \
  -4.25 4.25 0.0 1e-8 "$recurrence"
