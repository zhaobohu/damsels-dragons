#!/bin/bash
#SBATCH -D /path/to/project/annotation/genes/Anax_junius
#SBATCH -J AnaxJuni.gene.annot
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=48
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
conda activate /path/to/project/conda/envs/rna

THREADS="$SLURM_CPUS_PER_TASK"
DATA="/path/to/project/data"

# download RNA-seq reads
ACC1="ERR15165722"
#fasterq-dump -e "$THREADS" --split-files -p -t "$TMPDIR" -O . ${ACC1}

# trim RNA-seq reads with fastp
#fastp \
#  -i ${ACC1}_1.fastq -o ${ACC1}_1.trimmed.fastq \
#  -I ${ACC1}_2.fastq -O ${ACC1}_2.trimmed.fastq \
#  -w "$THREADS" \
#  --cut_tail \
#  --cut_tail_window_size 3 \
#  --cut_tail_mean_quality 25 \
#  --detect_adapter_for_pe \
#  --length_required 36 \
#  --trim_poly_x \
#  -j "${ACC1}.fastp.json" \
#  -h "${ACC1}.fastp.html"

# STAR index of the soft-masked genome
#STAR \
#  --runThreadN $THREADS \
#  --runMode genomeGenerate \
#  --genomeDir ./star.index \
#  --genomeFastaFiles $DATA/genomes/Anax.junius.clean.masked.fna \
#  --genomeSAindexNbases 14

# map trimmed reads with STAR (two-pass) and index the BAM
#STAR \
#  --runThreadN $THREADS \
#  --genomeDir ./star.index \
#  --readFilesIn "${ACC1}_1.trimmed.fastq" "${ACC1}_2.trimmed.fastq" \
#  --twopassMode Basic \
#  --outSAMtype BAM SortedByCoordinate \
#  --outSAMattributes All \
#  --outSAMstrandField intronMotif
#samtools index -@ $THREADS Aligned.sortedByCoord.out.bam

# writable copy of the AUGUSTUS config from the BRAKER3 container
#AUGUSTUS_CONFIG="$PWD/augustus_config"
#if [[ ! -d "$AUGUSTUS_CONFIG" ]]; then
#  singularity exec \
#    -B $PWD:$PWD \
#    /path/to/project/programs/braker3.sif \
#    cp -r /opt/Augustus/config "$AUGUSTUS_CONFIG"
#fi

# BRAKER3 run 1: RNA-seq BAM + OrthoDB v12 Arthropoda proteins
#singularity exec \
#  -B /path/to/project/annotation/genes/Anax_junius:/path/to/project/annotation/genes/Anax_junius \
#  -B /path/to/project:/path/to/project \
#  -B ${AUGUSTUS_CONFIG}:/opt/Augustus/config \
#  -B $HOME/.gm_key:/root/.gm_key \
#  /path/to/project/programs/braker3.sif braker.pl \
#    --genome=$DATA/genomes/Anax.junius.clean.masked.fna \
#    --bam=/path/to/project/annotation/genes/Anax_junius/Aligned.sortedByCoord.out.bam \
#    --prot_seq=$DATA/proteins/Arthropoda.OrthoDB.V12.faa \
#    --AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
#    --species=Anax.junius \
#    --threads $THREADS \
#    --workingdir=/path/to/project/annotation/genes/Anax_junius

# convert the BRAKER GTF to GFF3
#singularity exec \
#  -B /path/to/project:/path/to/project \
#  /path/to/project/programs/braker3.sif bash -c '/opt/Augustus/scripts/gtf2gff.pl --gff3 \
#    --out=/path/to/project/data/genes/Anax.junius.gff \
#    </path/to/project/annotation/genes/Anax_junius/braker.gtf'

# count genes on OZ* sequences; keep the longest isoform per gene
#awk '$3=="gene" && $1~/^OZ/' $DATA/genes/Anax.junius.gff \
#  | awk '{print $1}' | sort | uniq -c
#agat_sp_keep_longest_isoform.pl --gff $DATA/genes/Anax.junius.gff -o $DATA/genes/Anax.junius.longest.gff

# BUSCO on run-1 proteins (insecta_odb10)
#gffread -y braker.proteins.fna -g $DATA/genomes/Anax.junius.clean.masked.fna $DATA/genes/Anax.junius.gff
#conda activate /path/to/project/conda/envs/busco
#busco -i braker.proteins.fna -l insecta_odb10 -o busco -m proteins -c $THREADS
#conda activate /path/to/project/conda/envs/rna

# run-2 BRAKER evidence: OrthoDB + odonate proteins (species-prefixed), non-letters -> X
#pfx() { awk -v p="$1" '/^>/{sub(/^>/, ">" p "_")} 1' "$2"; }
#{
#  cat "$DATA/proteins/Arthropoda.OrthoDB.V12.faa"
#  pfx AeshJunc ../Aeshna_juncea/braker.aa
#  pfx ArchGran ../Archilestes_grandis/braker.aa
#  pfx HetaAmer ../../protein_evidence/Hetaerina.americana.proteins.faa
#  pfx IschEleg ../../protein_evidence/Ischnura.elegans.proteins.faa
#  pfx OrthCoer ../Orthetrum_coerulescens/braker.aa
#  pfx PlatPenn ../Platycnemis_pennipes/braker.aa
#  pfx SomaArct ../Somatochlora_arctica/braker.aa
#  pfx SympStri ../Sympetrum_striolatum/braker.aa
#} | awk '!/^>/{gsub(/[^A-Za-z]/,"X")} 1' > cat.prot.braker.faa

# GALBA evidence: odonate proteins only (species-prefixed)
#{ pfx AeshJunc ../Aeshna_juncea/braker.aa
#  pfx ArchGran ../Archilestes_grandis/braker.aa
#  pfx HetaAmer ../../protein_evidence/Hetaerina.americana.proteins.faa
#  pfx IschEleg ../../protein_evidence/Ischnura.elegans.proteins.faa
#  pfx OrthCoer ../Orthetrum_coerulescens/braker.aa
#  pfx PlatPenn ../Platycnemis_pennipes/braker.aa
#  pfx SomaArct ../Somatochlora_arctica/braker.aa
#  pfx SympStri ../Sympetrum_striolatum/braker.aa
#} > cat.prot.galba.faa

# BRAKER3 run 2: clear run-1 outputs, fresh AUGUSTUS config, combined protein evidence
#rm -rf GeneMark-ETP agat_log_Anax.junius Augustus species augustus_config
#AUGUSTUS_CONFIG="$PWD/augustus_config"
#if [[ ! -d "$AUGUSTUS_CONFIG" ]]; then
#  singularity exec \
#    -B $PWD:$PWD \
#    /path/to/project/programs/braker3.sif \
#    cp -r /opt/Augustus/config "$AUGUSTUS_CONFIG"
#fi
#singularity exec \
#  -B /path/to/project/annotation/genes/Anax_junius:/path/to/project/annotation/genes/Anax_junius \
#  -B /path/to/project:/path/to/project \
#  -B ${AUGUSTUS_CONFIG}:/opt/Augustus/config \
#  -B $HOME/.gm_key:/root/.gm_key \
#  /path/to/project/programs/braker3.sif braker.pl \
#    --genome=$DATA/genomes/Anax.junius.clean.masked.fna \
#    --bam=/path/to/project/annotation/genes/Anax_junius/Aligned.sortedByCoord.out.bam \
#    --prot_seq=/path/to/project/annotation/genes/Anax_junius/cat.prot.braker.faa \
#    --AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
#    --species=Anax.junius.re \
#    --threads $THREADS \
#    --workingdir=/path/to/project/annotation/genes/Anax_junius

# GALBA with odonate protein evidence only
#AUGUSTUS_CONFIG="$PWD/augustus_config"
#if [[ ! -d "$AUGUSTUS_CONFIG" ]]; then
#  singularity exec -B $PWD:$PWD \
#    /path/to/project/programs/galba.sif cp -r /opt/Augustus/config "$AUGUSTUS_CONFIG"
#fi
#rm -rf galba "$AUGUSTUS_CONFIG/species/Anax.junius.galba"
# scrub non-letter characters (. * etc.) from the evidence so DIAMOND accepts it
#awk '!/^>/{gsub(/[^A-Za-z]/,"X")} 1' cat.prot.galba.faa > tmp && mv tmp cat.prot.galba.faa
#singularity exec \
#  -B $PWD:$PWD -B /path/to/project:/path/to/project \
#  -B ${AUGUSTUS_CONFIG}:/opt/Augustus/config \
#  /path/to/project/programs/galba.sif galba.pl \
#    --species=Anax.junius.galba \
#    --genome=$DATA/genomes/Anax.junius.clean.masked.fna \
#    --prot_seq=$PWD/cat.prot.galba.faa \
#    --AUGUSTUS_CONFIG_PATH=/opt/Augustus/config \
#    --threads=$THREADS \
#    --workingdir=$PWD/galba

# re-run only GALBA's final DIAMOND filter step, on scrubbed evidence
#cd "$PWD/galba"
#awk '!/^>/{gsub(/[^A-Za-z]/,"X")} 1' prot_seq_all.fa > prot_seq_all.clean.fa
#singularity exec -B /path/to/project:/path/to/project \
#  /path/to/project/programs/galba.sif \
#  /opt/conda/bin/python3 /opt/GALBA/scripts/filter_gtf_by_diamond_against_ref.py \
#    -r prot_seq_all.clean.fa -g augustus.hints.gtf -o galba.gtf \
#    -a augustus.hints.aa -t 48 -d /opt/diamond
#cd ../

# merge BRAKER and GALBA gene sets with TSEBRA
#TSEBRA=/path/to/project/programs/TSEBRA
#python $TSEBRA/bin/rename_gtf.py --gtf braker.gtf      --prefix br --out br.gtf
#python $TSEBRA/bin/rename_gtf.py --gtf galba/galba.gtf --prefix ga --out ga.gtf
#python $TSEBRA/bin/tsebra.py \
#  -g br.gtf,ga.gtf \
#  -e hintsfile.gff,galba/hintsfile.gff \
#  -c $TSEBRA/config/default.cfg \
#  -o anax.combined.gtf

# rename merged genes, convert to GFF3, count genes, keep longest isoforms
#python /path/to/project/programs/TSEBRA/bin/rename_gtf.py \
#  --gtf anax.combined.gtf --prefix anax --out anax.renamed.gtf
#singularity exec \
#  -B /path/to/project:/path/to/project \
#  /path/to/project/programs/braker3.sif bash -c '/opt/Augustus/scripts/gtf2gff.pl --gff3 \
#    --out=/path/to/project/data/genes/Anax.junius.re.gff \
#    </path/to/project/annotation/genes/Anax_junius/anax.renamed.gtf'
#awk '$3=="gene" && $1~/^OZ/' $DATA/genes/Anax.junius.gff \
#  | awk '{print $1}' | sort | uniq -c
#agat_sp_keep_longest_isoform.pl --gff $DATA/genes/Anax.junius.gff -o $DATA/genes/Anax.junius.longest.gff

# BUSCO on longest-isoform proteins (insecta_odb12.2)
rm -rf busco_combined busco_downloads
agat_sp_extract_sequences.pl -g $DATA/genes/Anax.junius.longest.gff -f $DATA/genomes/Anax.junius.clean.fna -t cds -p -o anax.longest.aa
sed 's/\*$//' anax.longest.aa > anax.busco.aa
conda activate /path/to/project/conda/envs/busco
busco -i anax.busco.aa -l insecta_odb12.2 -o busco_combined -m proteins -c $THREADS
