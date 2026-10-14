#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.tree
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-10%10
# not /dev/null: failures before the script body are only reported here
#SBATCH -o slurm-%x.%A_%a.out
#SBATCH -e slurm-%x.%A_%a.err

# ten IQ-TREE searches from different seeds; 08b_support.sbatch picks the best by likelihood

set -eoux pipefail
mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
mkdir -p tmp
i="$SLURM_ARRAY_TASK_ID"
export TMPDIR="$(mktemp -d -p "$PWD/tmp" "iqtree.${SLURM_JOB_ID}.${i}.XXXX")"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/coitree
set -ux

THREADS="$SLURM_CPUS_PER_TASK"
outdir="run.${i}"; mkdir -p "$outdir"

# tree search for seed i; -g constrains the topology to the published backbone
iqtree \
  -p tree/partitions.nex \
  -g tree/constraint.tre \
  -m MFP+MERGE -nt "$THREADS" \
  -ninit 2000 -ntop 200 -nbest 20 -nstop 500 -pers 0.2 \
  -me 0.001 --seed $((0+i)) \
  -pre "${outdir}/run.${i}"
