#!/usr/bin/env bash
# P1.5: interval-width ladder at fixed maxdim=256 for the clean H0 + S.
# The P1.3 [-4.25, 4.25] run is omitted because it is already in progress.
# All intervals are validated against the dense initial spectrum
# [-3.983143..., 3.983143...]; they are diagnostics, not SCF bounds.
#SBATCH --job-name=qtt_aa_2d_p15
#SBATCH --output=p1-5-interval-%A_%a.out
#SBATCH --error=p1-5-interval-%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=12:00:00
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

HALF_WIDTHS=(4.05 5.0 6.0)
half_width="${HALF_WIDTHS[$((SLURM_ARRAY_TASK_ID - 1))]}"
label="$half_width"
label="${label//./p}"
output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside5_seed0p1_interval_ladder/task_0001/bound_${label}"

cd "$REPOSITORY_ROOT"
julia --startup-file=no --project=. \
  runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  1 256 "$output" "-$half_width" "$half_width"
