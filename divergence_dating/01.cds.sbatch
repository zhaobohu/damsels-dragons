#!/bin/bash
#SBATCH -D /path/to/project/divergence_dating
#SBATCH -J estimate
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
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
conda activate /path/to/project/conda/envs/pcm

THREADS="$SLURM_CPUS_PER_TASK"
DATA="/path/to/project/data"

# longest-isoform CDS of each species (AGAT)
for sp in Aeshna.juncea Anax.junius Archilestes.grandis Cordulegaster.boltonii \
  Ischnura.elegans Orthetrum.coerulescens Pantala.flavescens Platycnemis.pennipes \
  Somatochlora.arctica Sympetrum.striolatum; do
  [ -s ${sp}.longest.cds.fna ] && continue   # skip if already done
  /path/to/project/conda/envs/pcm/bin/agat_sp_extract_sequences.pl \
    -g $DATA/genes/${sp}.longest.gff \
    -f $DATA/genomes/${sp}.clean.masked.fna \
    -t cds \
    -o ${sp}.longest.cds.fna
done