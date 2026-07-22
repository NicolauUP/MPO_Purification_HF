#!/usr/bin/env bash
# Controlled chi=512 refinement of only V2=0, U=1 on the 32x32 square.
# 16 CPUs request about 32 GiB on GPP for MPO intermediate-memory headroom;
# Julia itself deliberately remains single-threaded.
#SBATCH --job-name=qtt_aa_2d_l5_c512
#SBATCH --output=chi512-%j.out
#SBATCH --error=chi512-%j.err
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

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/execute_square_campaign.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_chi512_clean.jl \
  1
