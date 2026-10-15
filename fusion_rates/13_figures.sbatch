#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.figures
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH -o slurm-%x.%A.out
#SBATCH -e slurm-%x.%A.err
#
# 13_figures.R on every map, then get_numbers.R clades (reads the 15_timeslice.R tables),
# so submit it after 15_timeslice.sbatch:
#   jid=$(sbatch --parsable 15_timeslice.sbatch | grep -oE '^[0-9]+' | tail -1)
#   sbatch --dependency=afterok:$jid 13_figures.sbatch
# (the grep drops the "Using default partition..." lines some clusters print with the job id)

set -eoux pipefail
mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"; echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/fusion
set -ux

# all maps -> fig_levels_permap.tsv; then the clade numbers for the text
Rscript 13_figures.R "$(Rscript -e 'cat(length(readRDS("simmaps.rds")))')"
Rscript get_numbers.R clades
