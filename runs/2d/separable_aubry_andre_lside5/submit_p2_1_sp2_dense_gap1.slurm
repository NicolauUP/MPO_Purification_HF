#!/usr/bin/env bash
# P2.1: compare fixed-Hamiltonian gap-1 SP2 with an independently constructed
# 1024x1024 exact occupied projector at maxdim 256 and 384.
#SBATCH --job-name=qtt_aa_2d_sp2_ed
#SBATCH --output=p2-1-sp2-ed-%A_%a.out
#SBATCH --error=p2-1-sp2-ed-%A_%a.err
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

MAXDIMS=(256 384)
maxdim="${MAXDIMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/square_lside5_sp2_dense_gap1/maxdim_${maxdim}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/compare_square_sp2_dense_projector.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  1 "$maxdim" "$output"
