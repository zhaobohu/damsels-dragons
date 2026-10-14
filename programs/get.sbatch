#!/bin/bash
#SBATCH -D /path/to/project/programs
#SBATCH -J get.sif.file
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -eoux pipefail

mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.err" >&2)

echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"
echo "Dir: $(pwd)"
echo "Start: $(date)"

rm -rf tmp
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +x; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/rna

# build the GALBA container image from Docker Hub
singularity build galba.sif docker://katharinahoff/galba-notebook:latest