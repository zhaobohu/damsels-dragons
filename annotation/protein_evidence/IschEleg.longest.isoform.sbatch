#!/bin/bash
#SBATCH -D /path/to/project/annotation/protein_evidence
#SBATCH -J IschEleg.longest.isoform
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=32
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

THREADS="$SLURM_CPUS_PER_TASK"
DATA="/path/to/project/data"

# count genes on OZ* sequences; keep the longest isoform per gene
#awk '$3=="gene" && $1~/^OZ/' $DATA/genes/Ischnura.elegans.gff \
#  | awk '{print $1}' | sort | uniq -c
#agat_sp_keep_longest_isoform.pl --gff $DATA/genes/Ischnura.elegans.gff -o $DATA/genes/Ischnura.elegans.longest.gff

# BUSCO on the full annotation's proteins (insecta_odb10)
#gffread -y Ischnura.elegans.cds.fna -g $DATA/genomes/Ischnura.elegans.clean.masked.fna $DATA/genes/Ischnura.elegans.gff
#conda activate /path/to/project/conda/envs/busco
#busco -i Ischnura.elegans.cds.fna -l insecta_odb10 -o busco -m proteins -c $THREADS
#conda activate /path/to/project/conda/envs/rna

# longest-isoform proteins, used as protein evidence for the odonate annotations
gffread -y Ischnura.elegans.proteins.faa -g ${DATA}/genomes/Ischnura.elegans.clean.fna ${DATA}/genes/Ischnura.elegans.longest.gff
