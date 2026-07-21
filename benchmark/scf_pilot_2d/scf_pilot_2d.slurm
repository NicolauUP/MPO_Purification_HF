#!/bin/bash
# Single 2D SCF pilot. Julia uses one thread; 16 CPUs request approximately
# 32 GiB on GPP, leaving headroom for the density, fields, and intermediates.
# Raise --cpus-per-task at submission only when a higher-memory cap requires it.
#SBATCH --job-name=qtt_scf_2d
#SBATCH --output=qtt-scf-2d-%j.out
#SBATCH --error=qtt-scf-2d-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail
module load julia/1.12.6

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0

BENCHMARK_ROOT="${MPO_BENCHMARK_ROOT:?set MPO_BENCHMARK_ROOT to an external result directory}"
OUTPUT_DIR="${BENCHMARK_ROOT}/scf_pilot_2d_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_DIR"

echo "Working directory: $WORKING_DIR"
echo "SCF pilot output: $OUTPUT_DIR"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK:-unknown} (Julia threads: $JULIA_NUM_THREADS)"
echo "Julia executable: $(command -v julia)"
julia --version

bash benchmark/scf_pilot_2d/run_scf_pilot_2d.sh "$OUTPUT_DIR" "$@"
