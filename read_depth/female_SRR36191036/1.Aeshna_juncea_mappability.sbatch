#!/bin/bash
#SBATCH -D /path/to/project/read_depth/female_SRR36191036
#SBATCH -J A.juncea.mappability
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Mappability mask (run once): windows.map.bed col4 = uniquely mappable fraction per window

set -eoux pipefail

mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)

echo "Start: $(date)"

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
WIN=5000        # window size (bp)
K=150           # k-mer length; raise toward read length to be less strict on long reads
E=1             # mismatches; >=1 flags divergent (not just identical) repeat copies
Y="OZ369555.1"  # juncea Y scaffold

# stop unless the Y scaffold is in the reference (else Y/microchromosome cross-mapping is invisible)
grep -q ">${Y}" "$REF" || { echo "ERROR: $Y absent from $REF -- add it before masking"; exit 1; }

# per-base mappability with GenMap (1 = uniquely mappable)
rm -rf gmidx
genmap index -F "$REF" -I gmidx
genmap map -K "$K" -E "$E" -I gmidx -O genmap -T "$THREADS" -bg

# uniquely mappable bases -> merged BED
awk '$4>=1' genmap.bedgraph | cut -f1-3 | sort -k1,1 -k2,2n \
  | bedtools merge -i - > mappable.bed

# fixed windows with the mappable fraction in col4
samtools faidx "$REF"
cut -f1,2 "${REF}.fai" > genome.txt
bedtools makewindows -g genome.txt -w "$WIN" > windows.bed
bedtools coverage -a windows.bed -b mappable.bed \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' > windows.map.bed

echo "windows: $(wc -l < windows.bed)  mappable intervals: $(wc -l < mappable.bed)"
echo "Done: $(date)"
