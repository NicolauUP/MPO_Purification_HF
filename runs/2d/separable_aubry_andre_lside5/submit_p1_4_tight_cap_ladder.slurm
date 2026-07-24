#!/usr/bin/env bash
# P1.4: cap ladder at the P1.3 tight fixed-Hamiltonian interval.
# Array task 1 uses maxdim=128; task 2 uses maxdim=384. The P1.3 maxdim=256
# control is intentionally not repeated here.
#SBATCH --job-name=qtt_aa_2d_p14
#SBATCH --output=p1-4-tight-cap-%A_%a.out
#SBATCH --error=p1-4-tight-cap-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=12:00:00
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

CAPS=(128 384)
maxdim="${CAPS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside5_seed0p1_tight_cap_ladder/task_0001/maxdim_${maxdim}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  1 "$maxdim" "$output" -4.25 4.25
