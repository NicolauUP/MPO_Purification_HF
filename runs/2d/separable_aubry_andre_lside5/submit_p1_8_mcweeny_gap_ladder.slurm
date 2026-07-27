#!/usr/bin/env bash
# P1.8: fixed-mu McWeeny gap ladder for the clean 32x32 square.
# Array tasks use gaps 0.5, 1.0, and 2.0, corresponding to checkerboard
# masses 0.25, 0.5, and 1.0. All numerical controls are otherwise identical.
#SBATCH --job-name=qtt_aa_2d_p18
#SBATCH --output=p1-8-gap-%A_%a.out
#SBATCH --error=p1-8-gap-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-3

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

TASK_LABELS=(gap_0p5 gap_1p0 gap_2p0)
task_label="${TASK_LABELS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/square_lside5_mcweeny_gap_ladder/task_$(printf '%04d' "$SLURM_ARRAY_TASK_ID")_${task_label}/maxdim_256"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/diagnose_square_mcweeny_mu_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_gap_ladder.jl \
  "$SLURM_ARRAY_TASK_ID" 256 "$output" \
  -4.25 4.25 0.0 1e-8
