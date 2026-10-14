#!/bin/bash
#SBATCH -D /path/to/project/palaeoptera_synteny/genespace_run
#SBATCH -J genespace
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
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
conda activate /path/to/project/conda/envs/genespace

export MCSCANX=/path/to/project/programs/MCScanX
THREADS="$SLURM_CPUS_PER_TASK"
DATA="/path/to/project/data"
GENOMES=/path/to/project/data/genomes
GENES=/path/to/project/data/genes/annevo

species=(
  Aeshna.juncea
  Somatochlora.meridionalis
  Somatochlora.arctica
  Serratella.ignita
  Anax.junius
  Sympetrum.sanguineum
  Tanypteryx.hageni
  Ischnura.elegans
  Cordulegaster.boltonii
  Pyrrhosoma.nymphula
  Orthetrum.coerulescens
  Anax.strenuus
  Sympetrum.striolatum
  Parafronurus.youi
  Ischnura.senegalensis
  Platycnemis.pennipes
  Ephemera.vulgata
  Ephemera.orientalis
  Ceriagrion.tenellum
  Ecdyonurus.torrentis
  Siphlonurus.alternatus
  Cloeon.dipterum
  Heteragrion.aequatoriale
  Baetis.rhodani
  Baetis.atlanticus
  Archilestes.grandis
  Procloeon.bifidum
  Hetaerina.titia
  Pantala.flavescens
)

# GENESPACE inputs: rawGenomeRepo/<id>/ with the longest-isoform gff and its proteins
# (id = species name with dots -> underscores)
#for s in "${species[@]}"; do
#  id="${s//./_}"
#  dir="rawGenomeRepo/${id}"
#  mkdir -p "$dir"
#  cp "${GENES}/${s}.longest.gff" "${dir}/${id}.gff"
#  gffread "${GENES}/${s}.longest.gff" \
#    -g "${GENOMES}/${s}.fna" \
#    -y "${dir}/${id}.protein.faa"
#done

# run GENESPACE
#Rscript genespace.R

# draw Fig. 1A
Rscript fig1a_riparian.R
