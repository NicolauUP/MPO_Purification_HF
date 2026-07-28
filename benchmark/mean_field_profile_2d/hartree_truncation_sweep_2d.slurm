#!/usr/bin/env bash
#SBATCH --job-name=qtt_hartree_sweep
#SBATCH --output=qtt-hartree-sweep-%j.out
#SBATCH --error=qtt-hartree-sweep-%j.err
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
OUTPUT="${ROOT}/hartree_truncation_sweep_lside10_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT"

/usr/bin/time -v -o "$OUTPUT/process_time.txt" \
  julia --startup-file=no --project=. \
  benchmark/mean_field_profile_2d/hartree_truncation_sweep_2d.jl \
  --output "$OUTPUT" \
  --side-level 10 \
  --preparation-cutoff 1e-10 \
  --preparation-maxdim 256 \
  --cutoffs 1e-12,1e-10,1e-8 \
  --maxdims 64,128,256,384,512 \
  --steps 50
