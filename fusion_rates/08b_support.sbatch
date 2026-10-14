#!/bin/bash
# no "#SBATCH -D": sbatch from inside the run directory
#SBATCH -J fusion.support
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
# not /dev/null: failures before the script body are only reported here
#SBATCH -o slurm-%x.%A_%a.out
#SBATCH -e slurm-%x.%A_%a.err

# UFBoot on the best 08_tree.sbatch topology. --wbtl is required: 09_date.R needs the
# bootstrap tree set. Check it holds N_BOOT trees before moving on:
#   grep -c ';' support/support.ufboot

set -eoux pipefail
mkdir -p logs support
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
rm -rf tmp; mkdir -p tmp; export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/coitree
set -ux

THREADS="$SLURM_CPUS_PER_TASK"

# scores of the ten searches; replicate count N_BOOT from config.R
grep -H "BEST SCORE" run.*/run.*.log | sort -t: -k3 -g
NBOOT="$(Rscript -e "source('config.R'); cat(N_BOOT)")"
: "${NBOOT:?could not read N_BOOT from config.R}"

# best-scoring run.N treefile, unless BEST_TREE is exported
BEST="${BEST_TREE:-$(grep -H "BEST SCORE" run.*/run.*.log | sort -t: -k3 -g | tail -1 | cut -d: -f1 | sed 's/\.log$/.treefile/')}"

# support on the fixed topology, writing the bootstrap trees (--wbtl)
iqtree -p tree/partitions.nex -g tree/constraint.tre -T "$THREADS" \
  -t "$BEST" \
  -B "$NBOOT" --wbtl -alrt "$NBOOT" --bnni \
  --seed 1 -pre support/support

grep -c ';' support/support.ufboot
