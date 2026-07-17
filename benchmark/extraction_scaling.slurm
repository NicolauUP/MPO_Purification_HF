#!/usr/bin/env bash
# Adjust CPU, memory, wall-time, module setup, partition, and account for the
# target cluster before submitting this file.
#SBATCH --job-name=qtt-extraction
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=qtt_extraction_%j.out
#SBATCH --error=qtt_extraction_%j.err

set -euo pipefail

# Example site-specific setup:
# module load julia

cd "${SLURM_SUBMIT_DIR:?submit this job from the MPO_Purification_HF directory}"
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

output_dir="benchmark_results/extraction_${SLURM_JOB_ID:-manual}"
bash benchmark/run_extraction_scaling.sh "$output_dir"
