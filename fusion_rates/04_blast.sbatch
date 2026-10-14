#!/bin/bash
# No "#SBATCH -D": submit from the run directory.
#SBATCH -J marker.blast
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH -o slurm-%x.%A_%a.out
#SBATCH -e slurm-%x.%A_%a.err

set -eoux pipefail
mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.err" >&2)
echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"
rm -rf tmp; mkdir -p tmp; export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/coitree
set -ux

THREADS="$SLURM_CPUS_PER_TASK"

# one reference per marker -> refs/<marker>.ref.fna
Rscript 04_find_refs.R

# BLAST each marker against its reference: screened COI records, or raw_<marker>.fasta
# (up to 20 HSPs per record; 05_select.R takes the union of their intervals)
mkdir -p blast
for ref in refs/*.ref.fna; do
  m="$(basename "$ref" .ref.fna)"
  if [ "$m" = "COI" ]; then q="screened.COI.fna"; else q="raw_${m}.fasta"; fi
  [ -s "$q" ] || { echo "skip $m (no query file)"; continue; }

  makeblastdb -in "$ref" -dbtype nucl -out "blast/${m}db" > /dev/null
  blastn -task blastn -query "$q" -db "blast/${m}db" \
         -max_target_seqs 1 -max_hsps 20 -num_threads "$THREADS" \
         -outfmt "6 qseqid qlen slen qstart qend sstart send pident" > "hits_${m}.tsv"
  echo "$m: $(wc -l < "hits_${m}.tsv") hits from $(grep -c '^>' "$q") queries"
done
