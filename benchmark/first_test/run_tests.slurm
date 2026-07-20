#!/usr/bin/env bash

# First cluster validation run for the maintained package test suite.
# Adjust account, QoS, resources, and module name only if the target cluster
# uses different values.
#SBATCH --job-name=qtt_hf_first_test
#SBATCH --output=qtt_hf_first_test-%j.out
#SBATCH --error=qtt_hf_first_test-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=112
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00

set -euo pipefail

module load julia/1.12.6

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit this job from the MPO_Purification_HF directory}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

echo "Working directory: $(pwd)"
echo "Julia executable: $(command -v julia)"
julia --version
echo "Julia threads: ${JULIA_NUM_THREADS}"

# Run the repository's maintained test entry point directly. Using Pkg.test()
# here would try to initialize/download the General registry on an offline
# compute node before the tests start.
julia --project=. test/runtests.jl
