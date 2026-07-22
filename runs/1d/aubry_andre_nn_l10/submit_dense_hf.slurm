#!/usr/bin/env bash
# Exact diagonalization of the same one-particle 1D HF equations as the MPO
# campaign. This is a three-task reference array, not many-body ED.
#SBATCH --job-name=dense_hf_aa_l10
#SBATCH --output=dense-hf-%A_%a.out
#SBATCH --error=dense-hf-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
  runs/common/execute_1d_dense_hf_campaign.jl \
  runs/1d/aubry_andre_nn_l10/campaign.jl \
  "$SLURM_ARRAY_TASK_ID"
