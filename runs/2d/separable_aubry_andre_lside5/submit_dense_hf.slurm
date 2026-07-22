#!/usr/bin/env bash
# Independent exact one-particle square HF reference for the same three cases.
#SBATCH --job-name=dense_hf_aa_2d_l5
#SBATCH --output=dense-hf-%A_%a.out
#SBATCH --error=dense-hf-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00
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

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/execute_square_dense_hf_campaign.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  "$SLURM_ARRAY_TASK_ID"
