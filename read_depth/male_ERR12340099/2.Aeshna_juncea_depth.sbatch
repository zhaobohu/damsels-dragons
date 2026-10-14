#!/bin/bash
#SBATCH -D /path/to/project/read_depth/male_ERR12340099
#SBATCH -J A.juncea.depth
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Windowed depth of uniquely mapped HiFi reads for one accession (male and female compared in R)

set -eoux pipefail

mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}.err" >&2)

echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"; echo "Dir: $(pwd)"; echo "Start: $(date)"

rm -rf tmp; mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

module load mamba/1.5.8
{ set +x; } 2>/dev/null
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/genome

THREADS=32
DATA="/path/to/project/data"
REF="${DATA}/genomes/Aeshna.juncea.clean.fna"
MAP="/path/to/project/read_depth/male_ERR12340099/windows.map.bed"  # from the mappability job
ACC1="ERR12340099"

# download and trim reads
#fasterq-dump -e "$THREADS" --split-files -p -t "$TMPDIR" -O . ${ACC1}
#fastp \
#  -i ${ACC1}.fastq -o ${ACC1}.trimmed.fastq -w "$THREADS" \
#  --cut_tail --cut_tail_window_size 1000 --cut_tail_mean_quality 30 \
#  --detect_adapter_for_pe --length_required 2000 --trim_poly_x \
#  -j "${ACC1}.fastp.json" -h "${ACC1}.fastp.html"

# map with minimap2; -p 0.1 -N 10 reports near-miss secondaries so ambiguous reads can be caught
minimap2 -ax map-hifi -p 0.1 -N 10 -t $(( THREADS - 6 )) \
  "$REF" ${ACC1}.trimmed.fastq | samtools sort -@ 6 -o align.bam -
samtools index align.bam

# names of reads with any secondary (256) or supplementary (2048) record
samtools view -@ "$THREADS" -f 256  align.bam | cut -f1 > tmp/m1
samtools view -@ "$THREADS" -f 2048 align.bam | cut -f1 > tmp/m2
sort -u tmp/m1 tmp/m2 > multi.qname

# primary records only, dropping every read in multi.qname
samtools view -@ "$THREADS" -h -F 0x900 -N ^multi.qname align.bam -o align.uniq.bam
samtools index align.uniq.bam

# mean depth per window (MAPQ >= 50); regions.bed.gz columns: chr start end frac depth
mosdepth -t "$THREADS" -Q 50 --by "$MAP" ${ACC1}.uniq align.uniq.bam

echo "multimapped reads removed: $(wc -l < multi.qname)"
echo "Done: $(date)"
