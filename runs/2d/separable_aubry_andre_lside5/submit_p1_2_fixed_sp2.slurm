#!/usr/bin/env bash
# P1.2: fixed-Hamiltonian SP2 cap sweep. No SCF or mean-field extraction.
# Tasks 1--3 correspond to the three campaign Hamiltonians. Each task runs
# maxdim 128, 256, 384, and 512 sequentially for a directly comparable trace.
# Eight allocated CPU cores provide the GPP default 16 GiB memory headroom;
# Julia remains single-threaded because this MPO path is not thread-parallel.
#SBATCH --job-name=qtt_aa_2d_p12
#SBATCH --output=p1-2-fixed-sp2-%A_%a.out
#SBATCH --error=p1-2-fixed-sp2-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
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

CAMPAIGN="runs/2d/separable_aubry_andre_lside5/campaign.jl"
TASK_INDEX="$SLURM_ARRAY_TASK_ID"
CAPS=(128 256 384 512)
RESULT_ROOT="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside5_seed0p1_fixed_sp2"

cd "$REPOSITORY_ROOT"
for maxdim in "${CAPS[@]}"; do
  output="${RESULT_ROOT}/task_$(printf '%04d' "$TASK_INDEX")/maxdim_${maxdim}"
  julia --startup-file=no --project=. \
    runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
    "$CAMPAIGN" "$TASK_INDEX" "$maxdim" "$output"
done
