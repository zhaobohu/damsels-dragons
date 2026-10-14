#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.date
#SBATCH --time=0-08:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
# not /dev/null: a launch failure is reported only here
#SBATCH -o slurm-%x.%A.out
#SBATCH -e slurm-%x.%A.err

# chronos on every bootstrap tree (not for a login node); 09_date.R resumes from its checkpoint

set -eoux pipefail
mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
mkdir -p tmp; export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/fusion
set -ux

Rscript 09_date.R
