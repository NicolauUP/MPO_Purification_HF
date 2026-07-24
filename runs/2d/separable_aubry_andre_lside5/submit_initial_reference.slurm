#!/usr/bin/env bash
# Exact dense spectra/projectors of H0 + S, before SP2 or SCF.
#SBATCH --job-name=qtt_aa_2d_initial
#SBATCH --output=initial-reference-%A_%a.out
#SBATCH --error=initial-reference-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=01:00:00
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
  runs/common/initial_square_hamiltonian_reference.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  "$SLURM_ARRAY_TASK_ID"
