#!/bin/bash
#SBATCH --job-name=qtt_hf_benchmark
#SBATCH --output=qtt_hf-%j.out
#SBATCH --error=qtt_hf-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail

module load julia

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit this job from the MPO_Purification_HF directory}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

RUN_NAME="extraction_${SLURM_JOB_ID}"
OUTPUT_DIR="${WORKING_DIR}/benchmark_results/${RUN_NAME}"

mkdir -p "$OUTPUT_DIR"

bash benchmark/run_extraction_scaling.sh "$OUTPUT_DIR"
