#!/bin/bash
#SBATCH --job-name=qtt_hf_2d_diag
#SBATCH --output=qtt-hf-2d-diag-%j.out
#SBATCH --error=qtt-hf-2d-diag-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=gp_ehpc
#SBATCH --account=epor78
#SBATCH --time=04:00:00

# Run the large smooth-density binary-carry Hartree diagnosis serially. Each case
# changes exactly one numerical control relative to the baseline: this is an
# accuracy diagnosis, not a parallel performance benchmark.
set -euo pipefail

module load julia/1.12.6

WORKING_DIR="${SLURM_SUBMIT_DIR:?submit from MPO_Purification_HF}"
cd "$WORKING_DIR"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/gpfs/projects/epor78/cluster_depot}"
export JULIA_PKG_PRECOMPILE_AUTO=0

DIAGNOSTIC_ROOT="${MPO_BENCHMARK_ROOT:?set MPO_BENCHMARK_ROOT to an external result directory}"
OUTPUT_DIR="${DIAGNOSTIC_ROOT}/binary_carry_hartree_diagnosis_${SLURM_JOB_ID}"
mkdir -p "$OUTPUT_DIR"

echo "Working directory: $WORKING_DIR"
echo "Diagnostic output: $OUTPUT_DIR"
echo "Julia executable: $(command -v julia)"
julia --version

run_case() {
  local label=$1
  shift
  echo "===== ${label} =====" | tee "$OUTPUT_DIR/${label}.stdout.txt"
  julia --startup-file=no --project=. \
    benchmark/extraction_scaling_2d/diagnose_tensorial_hartree.jl \
    --side-level 10 \
    --output "$OUTPUT_DIR/${label}.txt" \
    "$@" 2>&1 | tee -a "$OUTPUT_DIR/${label}.stdout.txt"
}

# Reference calculation: the settings used by the completed 2D benchmark.
run_case baseline --cutoff 1e-12 --maxdim 128

# Cutoff isolation: fixed maxdim=128.
run_case cutoff_1e-10 --cutoff 1e-10 --maxdim 128
run_case cutoff_1e-14 --cutoff 1e-14 --maxdim 128

# Bond-dimension isolation: fixed cutoff=1e-12.
run_case maxdim_256 --cutoff 1e-12 --maxdim 256
run_case maxdim_512 --cutoff 1e-12 --maxdim 512

echo "All binary-carry Hartree isolation cases completed."
