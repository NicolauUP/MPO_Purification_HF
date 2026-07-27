#!/usr/bin/env bash
# P2.3: compare intrinsic compression of the exact gap-2 projector with
# post-compression of independently generated SP2 projectors.
#SBATCH --job-name=qtt_aa_2d_projector_compress
#SBATCH --output=p2-3-projector-compress-%A_%a.out
#SBATCH --error=p2-3-projector-compress-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
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

MODES=(exact sp2 sp2)
MAXDIMS=(1024 256 384)
mode="${MODES[$((SLURM_ARRAY_TASK_ID - 1))]}"
maxdim="${MAXDIMS[$((SLURM_ARRAY_TASK_ID - 1))]}"
output="${MPO_RESULTS_ROOT}/square_lside5_projector_compression_gap2/${mode}_maxdim_${maxdim}"
cutoffs="1e-14,1e-12,1e-10,1e-8,1e-6"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/compare_square_projector_compression.jl \
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap1_gap2.jl \
  2 "$mode" "$maxdim" "$output" "$cutoffs"
