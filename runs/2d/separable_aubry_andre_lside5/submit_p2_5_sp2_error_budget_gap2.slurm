#!/usr/bin/env bash
# P2.5: fixed-Hamiltonian calibration of the proposed production error budget:
# operational cutoff 1e-10, idempotency 2e-4, relative trace error 1e-6.
#SBATCH --job-name=qtt_aa_2d_sp2_budget
#SBATCH --output=p2-5-sp2-budget-%A_%a.out
#SBATCH --error=p2-5-sp2-budget-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-2

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

MAXDIMS=(128 256)
maxdim="${MAXDIMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/square_lside5_sp2_error_budget_gap2/maxdim_${maxdim}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/compare_square_projector_compression.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  2 sp2 "$maxdim" "$output" \
  1e-8 1e-10 2e-4 1e-6
