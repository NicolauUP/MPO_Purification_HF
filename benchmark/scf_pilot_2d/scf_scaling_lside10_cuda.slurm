#!/usr/bin/env bash
#SBATCH --job-name=qtt_scf_l10
#SBATCH --output=qtt-scf-l10-cuda-%j.out
#SBATCH --error=qtt-scf-l10-cuda-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=acc
#SBATCH --gres=gpu:1
#SBATCH --qos=acc_ehpc
#SBATCH --account=epor78
#SBATCH --time=24:00:00

set -euo pipefail
bash benchmark/scf_pilot_2d/run_scf_scaling_cuda.sh 10
