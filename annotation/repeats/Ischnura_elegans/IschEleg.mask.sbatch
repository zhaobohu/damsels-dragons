#!/bin/bash
#SBATCH -D /path/to/project/annotation/repeats/Ischnura_elegans
#SBATCH -J IschEleg.mask
#SBATCH --time=3-00:00:00
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
/path/to/project/conda/envs/genome/share/RepeatMasker/famdb.py \
  -i $DATA/famdb \
  families \
  --format fasta_name \
  --include-class-in-name \
  -ad --curated 6961 \
  > Aeshnidae.curated.fa

# keep only the sequence ID in FASTA headers (once)
if [ ! -f $DATA/genomes/Ischnura.elegans.clean.fna ]; then
  awk '/^>/{print $1; next} {print}' $DATA/genomes/Ischnura.elegans.fna > $DATA/genomes/Ischnura.elegans.clean.fna
fi

# de novo repeat library with RepeatModeler (with LTR structure search)
BuildDatabase -name "Ischnura.elegans.wgs" $DATA/genomes/Ischnura.elegans.clean.fna
RepeatModeler -database "Ischnura.elegans.wgs" -threads $THREADS -LTRStruct

# combine de novo and curated libraries, dropping blank lines
cat Ischnura.elegans.wgs-families.fa Aeshnidae.curated.fa > temp.fna
sed '/^[[:space:]]*$/d' temp.fna > Ischnura.elegans.wgs.repeats.fna

# index the genome and soft-mask it with RepeatMasker
samtools faidx $DATA/genomes/Ischnura.elegans.clean.fna
RepeatMasker -xsmall -pa $THREADS -a -lib Ischnura.elegans.wgs.repeats.fna $DATA/genomes/Ischnura.elegans.clean.fna
