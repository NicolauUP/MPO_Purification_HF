#!/bin/bash
#SBATCH --job-name=qtt_hf_benchmark
#SBATCH --output=qtt_hf-%j.out
#SBATCH --error=qtt_hf-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
# The benchmark runs one sequential Julia process. Increase this only after
# establishing that a particular extraction kernel benefits from threading.
#SBATCH --cpus-per-task=1
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail

module load julia/1.12.6

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit this job from the MPO_Purification_HF directory}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0

RUN_NAME="extraction_${SLURM_JOB_ID}"
BENCHMARK_ROOT="${MPO_BENCHMARK_ROOT:?set MPO_BENCHMARK_ROOT to an external benchmark-result directory}"
OUTPUT_DIR="${BENCHMARK_ROOT}/${RUN_NAME}"

mkdir -p "$OUTPUT_DIR"

echo "Working directory: $WORKING_DIR"
echo "Benchmark output: $OUTPUT_DIR"
echo "Julia executable: $(command -v julia)"
julia --version
echo "Julia threads: $JULIA_NUM_THREADS"

bash benchmark/run_extraction_scaling.sh "$OUTPUT_DIR"
