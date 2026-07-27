#!/usr/bin/env bash
# P1.11: fixed-Hamiltonian canonical SP2 at checkerboard mass ±2 (gap Δ=4).
# No SCF or Hartree/Fock extraction is performed. The clean hopping spectrum
# and checkerboard anticommutation give |E| <= sqrt(20) < 4.75.
#SBATCH --job-name=qtt_aa_2d_sp2_gap4
#SBATCH --output=p1-11-sp2-gap4-%j.out
#SBATCH --error=p1-11-sp2-gap4-%j.err
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
  runs/2d/separable_aubry_andre_lside5/campaign_sp2_gap4.jl \
  1 256 \
  "${MPO_RESULTS_ROOT}/square_lside5_sp2_gap4/maxdim_256"
