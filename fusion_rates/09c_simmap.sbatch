#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.simmap
#SBATCH -p single
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH -o slurm-%x.%A_%a.out
#SBATCH -e slurm-%x.%A_%a.err
#
# One stochastic map per array task (09c_simmap_one.R). Submit as:
#   Rscript 09b_prep_simmap.R
#   Rscript 09d_make_todo.R <n_total>          # n_total is printed by 09b
#   sbatch --array=1-$(wc -l < todo.txt)%200 09c_simmap.sbatch
# Tasks past the end of todo.txt exit at once; %200 caps concurrency.

set -eoux pipefail
mkdir -p logs
echo "Job: ${SLURM_JOB_NAME} (${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

# the sampler is single-threaded; stop BLAS from grabbing cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/fusion
set -ux

Rscript 09c_simmap_one.R
