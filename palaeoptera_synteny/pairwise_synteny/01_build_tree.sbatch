#!/bin/bash
# Step 1 (MAFFT + IQ-TREE species tree) as a batch job: sbatch 01_build_tree.sbatch
#SBATCH -J figS1_tree
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem-per-cpu=2G
#SBATCH -t 48:00:00
#SBATCH -o figS1_tree_%j.out
#SBATCH -e figS1_tree_%j.err

set -euo pipefail

WORKDIR=/path/to/project/palaeoptera_synteny/pairwise_synteny
ENV_BIN=/path/to/project/conda/envs/phylo/bin

# put the conda env's bin on PATH (no conda hook in a batch shell)
if [ ! -d "$ENV_BIN" ]; then
  echo "ENV_BIN does not exist: $ENV_BIN" >&2
  echo "  find it with:  conda env list" >&2
  exit 1
fi
export PATH="$ENV_BIN:$PATH"
cd "$WORKDIR"

# check the tools are found
for prog in mafft iqtree3; do
  command -v "$prog" >/dev/null || { echo "$prog not found on PATH after adding $ENV_BIN" >&2; exit 1; }
done

# record job settings and tool paths
THREADS=${SLURM_NTASKS:-${SLURM_CPUS_ON_NODE:-8}}
echo "[job] $(date)  host $(hostname)  threads $THREADS"
echo "[job] mafft   $(command -v mafft)"
echo "[job] iqtree3 $(command -v iqtree3)"
echo "[job] python3 $(command -v python3)  $(python3 --version 2>&1)"

# --resume: a requeued job continues from the IQ-TREE checkpoint (starts fresh if tree/iqtree/ is empty)
python3 01_build_tree.py \
  --iqtree iqtree3 \
  --threads "$THREADS" \
  --resume

echo "[job] done $(date)"
ls -la tree/iqtree/tree.treefile
