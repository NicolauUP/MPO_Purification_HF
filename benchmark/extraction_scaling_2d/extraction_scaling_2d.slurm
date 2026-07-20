#!/bin/bash
#SBATCH --job-name=qtt_hf_2d_extract
#SBATCH --output=qtt-hf-2d-%j.out
#SBATCH --error=qtt-hf-2d-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
OUTPUT_DIR="${BENCHMARK_ROOT}/extraction_2d_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_DIR"

echo "Working directory: $WORKING_DIR"
echo "Benchmark output: $OUTPUT_DIR"
echo "Julia executable: $(command -v julia)"
julia --version

bash benchmark/extraction_scaling_2d/run_extraction_scaling_2d.sh "$OUTPUT_DIR" "$@"
