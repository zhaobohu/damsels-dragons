#!/bin/bash
#SBATCH -D /path/to/project/gene_fate/depth_check
#SBATCH -J sanity
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=64
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
conda activate /path/to/project/conda/envs/genome

# Version A: male HiFi reads -> hap1, per-gene CDS depth and heterozygosity
REF=/path/to/project/data/genomes/Aeshna.juncea.clean.masked.fna   # ioAesJunc2.hap1.1, all scaffolds
READS=ERR12340099.trimmed.fastq                                     # HiFi reads of the assembled male
GFF=/path/to/project/data/genes/Aeshna.juncea.longest.gff
T=64
MAPQ=20        # rerun with MAPQ=0 once to see how much depth is lost to gametolog ambiguity

# 1. map
samtools faidx $REF
minimap2 -t $T -ax map-hifi --secondary=no $REF $READS \
  | samtools sort -@ 8 -m 2G -o male.hap1.bam -
samtools index male.hap1.bam

# 2. CDS bed, one transcript per gene
agat_sp_keep_longest_isoform.pl --gff $GFF -o longest.gff3
# assumes BRAKER-style ids (gene g1, transcript g1.t1). Check `head cds.bed`.
awk -F'\t' '$3=="CDS"{ match($9,/Parent=[^;]+/); p=substr($9,RSTART+7,RLENGTH-7);
                       sub(/\.t[0-9]+$/,"",p); print $1"\t"$4-1"\t"$5"\t"p }' longest.gff3 \
  | sort -k1,1 -k2,2n > cds.bed
echo "genes: $(cut -f4 cds.bed | sort -u | wc -l)"

# 3. depth
mosdepth -t 4 -n --mapq $MAPQ --by cds.bed male.cds male.hap1.bam   # per-CDS-interval mean
mosdepth -t 4 -n --mapq $MAPQ --by 200000 male.win male.hap1.bam    # 200 kb windows (median depth, step 4)

# 4. callable CDS bp per gene (depth within 0.25x-2x of the genome median)
MED=$(zcat male.win.regions.bed.gz | cut -f4 | sort -g | awk '{a[NR]=$1} END{print a[int(NR/2)]}')
LO=$(awk -v m=$MED 'BEGIN{printf "%d", m*0.25}')
HI=$(awk -v m=$MED 'BEGIN{printf "%d", m*2}')
echo "window median depth $MED; callable range $LO-$HI"
mosdepth -t 4 -n --mapq $MAPQ --quantize 0:$LO:$HI: male.q male.hap1.bam
zcat male.q.quantized.bed.gz | awk -v r="$LO:$HI" '$4==r' > callable.bed
bedtools intersect -a cds.bed -b callable.bed \
  | awk '{n[$4]+=$3-$2} END{for(g in n) print g"\t"n[g]}' > callable_cds_bp.tsv

# 5. heterozygous SNPs in CDS
bcftools mpileup -f $REF -X pacbio-ccs -q $MAPQ -a FORMAT/AD,FORMAT/DP \
    -R cds.bed --threads 8 male.hap1.bam \
  | bcftools call -mv -Oz -o male.cds.vcf.gz
bcftools index male.cds.vcf.gz
bcftools view -v snps -i 'GT="het" && QUAL>=20 && FMT/DP>=5' male.cds.vcf.gz \
  | bcftools query -f '%CHROM\t%POS0\t%POS\t[%AD]\n' > het_sites.bed
bedtools intersect -a het_sites.bed -b cds.bed -wb \
  | awk '{n[$8]++} END{for(g in n) print g"\t"n[g]}' > het_per_gene.tsv

echo "done. now: Rscript 02_build_table.R"
