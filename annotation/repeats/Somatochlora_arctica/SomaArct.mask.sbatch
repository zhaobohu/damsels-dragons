#!/bin/bash
#SBATCH -D /path/to/project/annotation/repeats/Somatochlora_arctica
#SBATCH -J SomaArct.mask
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=32
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

module load mamba/1.5.8
{ set +x; } 2>/dev/null
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/repeatmodeler

THREADS="$SLURM_CPUS_PER_TASK"
DATA="/path/to/project/data"

# Dfam curated library: fetch FamDB partitions 0 and 15 once, export taxon 6961 families
if [ ! -f $DATA/famdb/dfam39_full.0.h5 ]; then
  wget -P $DATA/famdb https://www.dfam.org/releases/current/families/FamDB/dfam39_full.0.h5.gz
  wget -P $DATA/famdb https://www.dfam.org/releases/current/families/FamDB/dfam39_full.15.h5.gz
  gunzip $DATA/famdb/dfam39_full.0.h5.gz $DATA/famdb/dfam39_full.15.h5.gz
fi
#/path/to/project/conda/envs/genome/share/RepeatMasker/famdb.py \
#  -i $DATA/famdb \
#  families \
#  --format fasta_name \
#  --include-class-in-name \
#  -ad --curated 6961 \
#  > Aeshnidae.curated.fa

# keep only the sequence ID in FASTA headers (once)
if [ ! -f $DATA/genomes/Somatochlora.arctica.clean.fna ]; then
  awk '/^>/{print $1; next} {print}' $DATA/genomes/Somatochlora.arctica.fna > $DATA/genomes/Somatochlora.arctica.clean.fna
fi

# de novo repeat library with RepeatModeler (with LTR structure search)
# -recoverDir resumes the interrupted RepeatModeler run
#BuildDatabase -name "Somatochlora.arctica.wgs" $DATA/genomes/Somatochlora.arctica.clean.fna
RepeatModeler -database "Somatochlora.arctica.wgs" -threads $THREADS -LTRStruct -recoverDir RM_3860440.ThuMay281442582026

# combine de novo and curated libraries, dropping blank lines
cat Somatochlora.arctica.wgs-families.fa Aeshnidae.curated.fa > temp.fna
sed '/^[[:space:]]*$/d' temp.fna > Somatochlora.arctica.wgs.repeats.fna

# index the genome and soft-mask it with RepeatMasker
samtools faidx $DATA/genomes/Somatochlora.arctica.clean.fna
RepeatMasker -xsmall -pa $THREADS -a -lib Somatochlora.arctica.wgs.repeats.fna $DATA/genomes/Somatochlora.arctica.clean.fna
