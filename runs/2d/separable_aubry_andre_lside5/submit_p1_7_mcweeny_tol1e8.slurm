#!/usr/bin/env bash
# P1.7: direct-McWeeny cutoff control. This differs from P1.6 only by
# itensors_tol=1e-8 versus the campaign's 1e-14.
#SBATCH --job-name=qtt_aa_2d_p17
#SBATCH --output=p1-7-mcweeny-tol1e8-%j.out
#SBATCH --error=p1-7-mcweeny-tol1e8-%j.err
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
  runs/common/diagnose_square_mcweeny_mu_fixed_hamiltonian.jl \
  runs/2d/separable_aubry_andre_lside5/campaign.jl \
  1 256 \
  "${MPO_RESULTS_ROOT}/separable_aubry_andre_lside5_seed0p1_mcweeny_mu_tol1e8/task_0001/maxdim_256" \
  -4.25 4.25 0.0 1e-8
