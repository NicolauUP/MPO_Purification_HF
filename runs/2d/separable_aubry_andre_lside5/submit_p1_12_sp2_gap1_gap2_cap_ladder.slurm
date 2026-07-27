#!/usr/bin/env bash
# P1.12: fixed-Hamiltonian canonical SP2 at gaps 1 and 2, each at
# maxdim 256 and 384. No SCF or Hartree/Fock extraction is performed.
#SBATCH --job-name=qtt_aa_2d_sp2_midgap
#SBATCH --output=p1-12-sp2-midgap-%A_%a.out
#SBATCH --error=p1-12-sp2-midgap-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00
#SBATCH --array=1-4

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

CAMPAIGN_TASKS=(1 1 2 2)
GAP_LABELS=(gap_1p0 gap_1p0 gap_2p0 gap_2p0)
MAXDIMS=(256 384 256 384)
index=$((SLURM_ARRAY_TASK_ID - 1))
campaign_task="${CAMPAIGN_TASKS[$index]}"
gap_label="${GAP_LABELS[$index]}"
maxdim="${MAXDIMS[$index]}"
output="${MPO_RESULTS_ROOT}/square_lside5_sp2_gap1_gap2/${gap_label}/maxdim_${maxdim}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  "$campaign_task" "$maxdim" "$output"
