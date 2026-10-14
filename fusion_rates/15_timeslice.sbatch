#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.timeslice
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G
#SBATCH -o slurm-%x.%A.out
#SBATCH -e slurm-%x.%A.err
#
# 15_timeslice.R forks one worker per core; memory ~ one copy of simmaps.rds plus a little per worker

set -eoux pipefail
mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"; echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

# one thread per forked worker
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/fusion
set -ux

Rscript 15_timeslice.R
