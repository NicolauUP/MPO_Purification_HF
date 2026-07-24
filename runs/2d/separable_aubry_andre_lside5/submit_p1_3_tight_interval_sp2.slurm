#!/usr/bin/env bash
# P1.3: only the initial-SP2 spectral scaling interval is changed.
# The exact dense spectrum of task 1 is [-3.98315, 3.98315], so [-4.25, 4.25]
# safely encloses H0 + S with 0.2668 energy units of margin on either side.
# This is a fixed-Hamiltonian diagnostic, not an SCF production bound.
#SBATCH --job-name=qtt_aa_2d_p13
#SBATCH --output=p1-3-tight-sp2-%j.out
#SBATCH --error=p1-3-tight-sp2-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=08:00:00

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
  runs/common/diagnose_square_sp2_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  1 256 \
  "${MPO_RESULTS_ROOT}/separable_aubry_andre_lside5_seed0p1_tight_sp2/task_0001/maxdim_256" \
  -4.25 4.25
