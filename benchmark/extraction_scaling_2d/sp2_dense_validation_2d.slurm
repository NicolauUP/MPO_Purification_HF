#!/bin/bash
#SBATCH --job-name=qtt_sp2_dense_check
#SBATCH --output=qtt-sp2-dense-check-%j.out
#SBATCH --error=qtt-sp2-dense-check-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00

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
OUTPUT_DIR="${BENCHMARK_ROOT}/sp2_dense_validation_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_DIR"

echo "Working directory: $WORKING_DIR"
echo "Validation output: $OUTPUT_DIR"
echo "Julia executable: $(command -v julia)"
julia --version

bash benchmark/extraction_scaling_2d/run_sp2_dense_validation_2d.sh "$OUTPUT_DIR" "$@"
