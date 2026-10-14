#!/bin/bash
#SBATCH -J ncbi_dl
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -o ncbi_dl_%j.out
#SBATCH -e ncbi_dl_%j.err

set -euo pipefail
DATASETS=/path/to/project/palaeoptera_synteny/download_genomes/datasets
WORK=/path/to/project/palaeoptera_synteny/download_genomes
LEVELS=chromosome,complete
export NCBI_API_KEY="${NCBI_API_KEY:?export your NCBI API key before running}"
cd "$WORK"

# get the NCBI datasets CLI
#curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
#chmod +x datasets

# dehydrated download of all chromosome-level Palaeoptera (taxon 33339) GenBank genomes
#"$DATASETS" download genome taxon 33339 \
#  --assembly-level "$LEVELS" \
#  --assembly-source genbank \
#  --include genome \
#  --dehydrated \
#  --no-progressbar \
#  --filename palaeoptera.zip

#unzip -o palaeoptera.zip -d raw

# rehydrate (fetch the sequences), retrying up to 5 times
#for i in 1 2 3 4 5; do
#  "$DATASETS" rehydrate --directory raw --max-workers 6 --gzip \
#    --api-key "$NCBI_API_KEY" --no-progressbar && break
#  echo "attempt $i incomplete, sleeping 120s"
#  sleep 120
#done

# collect the genome fasta files in fasta/
#mkdir -p fasta
#find raw/ncbi_dataset/data -name '*_genomic.fna*' -exec mv {} fasta/ \;
#ls fasta | wc -l
#du -sh fasta

# decompress
cd /path/to/project/palaeoptera_synteny/download_genomes/fasta
gunzip *.fna.gz
