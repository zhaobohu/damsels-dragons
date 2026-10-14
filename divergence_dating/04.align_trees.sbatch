#!/bin/bash
#SBATCH -D /path/to/project/divergence_dating
#SBATCH -J align
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# PRANK codon alignment, then IQ-TREE gene trees, both seeded
set -eoux pipefail

SEED=67   # unseeded PRANK placed gaps differently between runs on identical input

mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.err" >&2)

echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"
echo "Dir: $(pwd)"
echo "Seed: $SEED"
echo "Start: $(date)"

rm -rf tmp
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +x; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/pcm

shopt -s nullglob

# codon alignment with PRANK, up to 3 tries per locus
for f in per_locus/*/*.fasta; do
  d=$(dirname $f); b=$(basename $f .fasta)
  [ -f $d/$b.best.fas ] && continue
  ok=0
  for try in 1 2 3; do
    if prank -d=$f -o=$d/$b -codon -F -seed=$SEED; then ok=1; break; fi
    echo "PRANK attempt $try failed: $b" >&2
    rm -f $d/$b.best.fas $d/$b.best.dnd
  done
  [ $ok -eq 1 ] || echo "PRANK GAVE UP: $b" >&2
done
echo "aligned: $(ls per_locus/*/*.best.fas | wc -l) / $(ls per_locus/*/*.fasta | wc -l)"

# ML gene trees with IQ-TREE, up to 3 tries per locus
for f in per_locus/*/*.best.fas; do
  d=$(dirname $f); b=$(basename $f .best.fas)
  [ -f $d/$b.treefile ] && continue
  ok=0
  for try in 1 2 3; do
    if iqtree3 -s $f -m MFP -B 1000 -T 1 -seed $SEED --prefix $d/$b -redo; then
      ok=1; break
    fi
    echo "IQ-TREE attempt $try failed: $b" >&2
    rm -f $d/$b.treefile $d/$b.ckp.gz
  done
  [ $ok -eq 1 ] || echo "IQ-TREE GAVE UP: $b" >&2
done
echo "trees: $(ls per_locus/*/*.treefile | wc -l) / $(ls per_locus/*/*.best.fas | wc -l)"

# per-locus alignment lengths and md5s, to compare against future runs
for f in per_locus/*/*.best.fas; do
  b=$(basename $f .best.fas)
  L=$(awk '/^>/{if(n){print n; exit}; next}{n+=length($0)}' $f)
  printf '%s\t%s\n' "$b" "$L"
done | tee alignment_lengths.tsv | awk -F'\t' '{s+=$2} END{print "total aligned nt:", s}'

md5sum per_locus/*/*.best.fas > alignment_md5.txt
echo "wrote alignment_lengths.tsv and alignment_md5.txt"
