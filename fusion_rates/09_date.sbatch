#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.date
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
# not /dev/null: a launch failure is reported only here
#SBATCH -o slurm-%x.%A.out
#SBATCH -e slurm-%x.%A.err

# LSD2 on the ML tree and every bootstrap tree (about 1 min per tree per core; not for a login node)
# LSD2 binary: set LSD2 below, or export LSD2=... before sbatch

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

export LSD2="${LSD2:-/path/to/project/conda/envs/lsd2/bin/lsd2}"
Rscript 09_date.R
