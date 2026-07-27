#!/usr/bin/env bash
# P2.2: determine whether gap-1 SP2 at maxdim 384 can improve beyond the
# default 1e-3 stopping threshold before MPO truncation causes stagnation.
#SBATCH --job-name=qtt_aa_2d_sp2_strict
#SBATCH --output=p2-2-sp2-strict-%j.out
#SBATCH --error=p2-2-sp2-strict-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

output="${MPO_RESULTS_ROOT}/square_lside5_sp2_dense_gap1_strict/maxdim_384_tol_1e-4"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/compare_square_sp2_dense_projector.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  1 384 "$output" 1e-4
