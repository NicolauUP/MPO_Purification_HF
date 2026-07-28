#!/usr/bin/env bash
#SBATCH --job-name=qtt_mf_profile
#SBATCH --output=qtt-mf-profile-%j.out
#SBATCH --error=qtt-mf-profile-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail
module load julia/1.12.6

cd "${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0

ROOT="${MPO_BENCHMARK_ROOT:?set MPO_BENCHMARK_ROOT to an external result directory}"
OUTPUT="${ROOT}/mean_field_component_profile_lside10_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT"

/usr/bin/time -v -o "$OUTPUT/process_time.txt" \
  julia --startup-file=no --project=. \
  benchmark/mean_field_profile_2d/component_profile_2d.jl \
  --output "$OUTPUT" \
  --side-level 10 \
  --cutoff 1e-10 \
  --maxdim 256 \
  --steps 50
