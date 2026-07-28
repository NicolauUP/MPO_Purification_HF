#!/usr/bin/env bash
# Diagonalize the N=4096 initial Hamiltonian once, then compress its exact
# occupied projector at a common cutoff and four QTT caps.
#SBATCH --job-name=qtt_aa_l6_exact_projector
#SBATCH --output=aa-l6-exact-projector-%j.out
#SBATCH --error=aa-l6-exact-projector-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail
module load julia/1.12.6

REPOSITORY_ROOT="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=32
export MKL_NUM_THREADS=32
export MPO_RESULTS_ROOT="${MPO_RESULTS_ROOT:?set an external result root before submitting}"

output="${MPO_RESULTS_ROOT}/separable_aubry_andre_lside6_exact_projector/cutoff_1e-8"
cd "$REPOSITORY_ROOT"
mkdir -p "$(dirname "$output")"
echo "host=$(hostname)"
echo "exact-projector output=$output"

/usr/bin/time -v -o "$(dirname "$output")/process_time.txt" \
  julia --startup-file=no --project=. \
  runs/common/compress_exact_square_projector_local.jl \
  runs/2d/separable_aubry_andre_lside6/campaign.jl \
  1 128,256,512,1024 1e-8 "$output" 32
