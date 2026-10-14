#!/usr/bin/env bash
#SBATCH -D /path/to/project/microchromosome/sex_kmers
#SBATCH -J kmer_sex
#SBATCH -p single
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=2G               # = 128G total; meryl-lookup holds the k-mer DB in RAM
#SBATCH -t 2-00:00:00
#SBATCH -o kmer_sex.%j.out
#SBATCH -e kmer_sex.%j.err
# Sex-specific 31-mer analysis of A. juncea on the male hap1 assembly.
# Submit from the directory that contains summarize_kmer_sex.R: sbatch kmer.sh
set -euo pipefail

{ set +ux; }
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/kmer
set -ux

# config
K=31
W=100000                              # window size (bp)
THREADS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-16}}   # SLURM_CPUS_PER_TASK was unset in job 1024750
MEM_GB=100                            # meryl count memory cap; keep BELOW the SLURM allocation (128G), not equal to it

# assemblies and reads (assembly identities checked against the NCBI sequence report in job 1024750)
HAP1=/path/to/project/data/genomes/Aeshna.juncea.fna                    # male hap1, chromosome-level (GCA_977065675.1)
HAP2=/path/to/project/data/genomes/Aeshna.juncea.2.1.fna                # male hap2, contig-level
FEM=/path/to/project/data/genomes/Aeshna.juncea.female.AK7.fna         # female primary (GCA_053455795.1), isolate AK7
MALE_READS=/path/to/project/data/reads/ERR12340099.fastq                # space-separated list is fine; leave unquoted where used
FEM_READS=/path/to/project/data/reads/SRR36191036.fastq
MM_PRESET_MALE=map-hifi               # minimap2 preset per sex: map-hifi / map-ont / sr
MM_PRESET_FEM=map-hifi

# chromosome selection (step 0): sequence ID, NCBI chromosome name (e.g. "14") or auto (= shortest molecule)
ACC_HAP1=GCA_977065675.1              # accession the HAP1 fasta came from; used only to fetch the sequence report
ACC_FEM=GCA_053455795.1
MICRO=OZ369553.1                      # chromosome 12, 48.2 Mb - the smallest non-sex chromosome
NEOX=OZ369554.1                       # SUPER_X, 184.5 Mb
NEOY=OZ369555.1                       # SUPER_Y, 166.7 Mb - male-limited, so it MUST NOT sit in the autosomal null
EXCLUDE="OZ369556.1"                  # dropped entirely (MT: 16 kb, wildly different copy number)
# the fused X is split into the neo-X arm (has a neo-Y counterpart) and the ancestral X arm
NEOX_BREAK=$(tr -dc '0-9' < /path/to/project/gene_pairs_ds/breakpoint.txt)   # fusion breakpoint on NEOX, bp (from the pairs run). Empty = treat the X as one block.
NEOX_ARM=left                         # which side of the breakpoint is the NEO-X arm: left | right
MIN_CHR_LEN=10000000                  # fallback only, used if step 0 metadata is unavailable
DATASETS=datasets                     # or /path/to/project/palaeoptera_synteny/download_genomes/datasets
export NCBI_API_KEY=${NCBI_API_KEY:-} # optional; set it if the summary call gets rate-limited

MIN_MALE=auto                         # solid-k-mer threshold for male reads; auto = first valley of the histogram
MIN_FEM=auto
N_HET_AUTO=3                          # how many (largest) autosomes to call alongside the candidates in step 9
# step 9 depth band: het calls and callable sites are both counted only within it (ceiling drops collapsed repeats)
HET_DP_MIN=8
HET_DP_MAX=40
NPERM=10000                           # block permutations in the R summary

OUT=kmer_sex_out
# steps: 0 chrom classes, 1 count, 3 set operations, 4 k-mer windows, 5 read depth, 6 assembly alignments,
#        7 length/GC, 8 R summary, 9 heterozygosity (can also be submitted alone)
STEPS="0 1 3 4 5 6 7 8 9"
# optional fetch:  datasets download genome accession GCA_977065675.1 --include genome --filename hap1.zip

# summarize_kmer_sex.R is looked for in the submit dir first, then the -D dir
SCRIPT_DIR=${SLURM_SUBMIT_DIR:-$(pwd)}
[[ -f "$SCRIPT_DIR/summarize_kmer_sex.R" ]] || SCRIPT_DIR=$(pwd)
run_step() { [[ " $STEPS " == *" $1 "* ]]; }

# distinct k-mers in a meryl DB (parses 'meryl statistics'; check the output format if this returns empty)
nk() { meryl statistics "$1" | awk '$1=="distinct"{print $2; exit}'; }

# meryl histogram -> "valley peak": first local minimum, then the highest bin after it
hist_valley_peak() {
  awk 'BEGIN{v=0}{ if(v==0 && NR>1 && $2>prev) v=pv; if(v>0 && $2>pk){pk=$2; p=$1}; prev=$2; pv=$1 } END{print v, p+0}' "$1"
}

# place several k-mer DBs on a fasta in one pass and count k-mer start positions per window
# (per-k-mer BED streamed through a FIFO; writes tables/windows_<label>.tsv for each label)
#   usage: lookup_windows <fasta> <db1> <label1> [<db2> <label2> ...]
lookup_windows() {
  local fasta=$1; shift
  local dbs=() labels=()
  while [[ $# -gt 0 ]]; do dbs+=("$1"); labels+=("$2"); shift 2; done
  # FIFO goes in local scratch, not on /work: named pipes are unreliable on parallel filesystems
  local fifo; fifo=$(mktemp -u -p "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}" lookup.XXXXXX.fifo); mkfifo "$fifo"
  ( awk -v W=$W -v PFX="tables/windows_" -v OFS="\t" \
        '{b=int($2/W); c[$4 SUBSEP $1 SUBSEP b]++}
         END{for(k in c){split(k,a,SUBSEP); print a[2], a[3]*W, (a[3]+1)*W, c[k] > (PFX a[1] ".raw")}}' < "$fifo" ) &
  local pid=$!
  meryl-lookup -bed -sequence "$fasta" -mers "${dbs[@]}" -labels "${labels[@]}" -output "$fifo"
  wait $pid
  rm -f "$fifo"
  local l
  for l in "${labels[@]}"; do
    [[ -s tables/windows_$l.raw ]] || { echo "ERROR: tables/windows_$l.raw is empty"; exit 1; }
    sort -k1,1 -k2,2n tables/windows_$l.raw > tables/windows_$l.tsv && rm -f tables/windows_$l.raw
  done
}

# fraction of each sequence covered by PAF alignments
paf_cov() {   # $1 paf  $2 P|ALL (primary only / all incl. secondary)  $3 target|query  $4 genome file (name len)  $5 out
  awk -v mode=$2 -v side=$3 -v OFS="\t" '{
      tp="?"; for(i=13;i<=NF;i++) if($i ~ /^tp:A:/) tp=substr($i,6,1)
      if(mode=="P" && tp!="P") next
      if(side=="target") print $6,$8,$9; else print $1,$3,$4 }' "$1" \
  | sort -k1,1 -k2,2n | bedtools merge -i - \
  | awk -v OFS="\t" 'NR==FNR{len[$1]=$2; next} {cov[$1]+=$3-$2}
         END{for(c in len) print c, len[c], cov[c]+0, (cov[c]+0)/len[c]}' "$4" - \
  | sort -k1,1 > "$5"
}

# PAF -> fixed 13 columns (12 standard + tp tag) for R
paf_slim() {
  awk -v OFS="\t" '{tp="?"; for(i=13;i<=NF;i++) if($i ~ /^tp:A:/) tp=substr($i,6,1);
                    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,tp}' "$1" > "$2"
}

# output folders
mkdir -p $OUT/{meryl,tables,beds,aln,reads}
cd $OUT
echo "K=$K W=$W THREADS=$THREADS MICRO=$MICRO NEOX=$NEOX STEPS=$STEPS"

# step 0: tables/chrom_classes.tsv from the NCBI sequence report; resolve MICRO / NEOX / NEOY to sequence IDs
if run_step 0; then
  $DATASETS summary genome accession $ACC_HAP1 --report sequence --as-json-lines > tables/hap1_seqreport.jsonl
  [[ -s tables/hap1_seqreport.jsonl ]] || { echo "ERROR: empty sequence report for $ACC_HAP1 (no network on the compute node? wrong accession?)"; exit 1; }
  echo "--- first sequence-report record (check the keys if the table below looks wrong):"
  head -1 tables/hap1_seqreport.jsonl
  # JSON keys read tolerantly: the sequence-report schema differs between datasets releases
  jq -r '[(.genbank_accession // .accession_version // .accession // .refseq_accession // ""),
          (.chr_name // .chromosome // .assigned_molecule // ""),
          (.role // .sequence_role // ""),
          (.length // .seq_length // .sequence_length // 0),
          (.sequence_name // "")] | @tsv' \
      tables/hap1_seqreport.jsonl > tables/hap1_seqmeta.tsv
  awk -F'\t' '$1==""||$4==0{bad++} END{if(bad>0){print "ERROR: "bad" records with no accession or length - inspect the JSON keys above" > "/dev/stderr"; exit 1}}' tables/hap1_seqmeta.tsv

  # do the metadata accessions actually match the FASTA headers?
  seqkit fx2tab -n -i -l $HAP1 > tables/hap1.genome
  n_match=$(awk -F'\t' 'NR==FNR{m[$1]; next} ($1 in m){n++} END{print n+0}' tables/hap1_seqmeta.tsv tables/hap1.genome)
  n_fa=$(wc -l < tables/hap1.genome)
  echo "sequence IDs matching metadata: $n_match / $n_fa"
  if ! (( n_match > 0 )); then
    echo "ERROR: no FASTA sequence ID matches the sequence report for $ACC_HAP1."
    echo "  HAP1 fasta        : $HAP1"
    echo "  first FASTA IDs   : $(cut -f1 tables/hap1.genome | head -3 | tr '\n' ' ')"
    echo "  first metadata IDs: $(cut -f1 tables/hap1_seqmeta.tsv | head -3 | tr '\n' ' ')"
    echo "  metadata seq names: $(cut -f5 tables/hap1_seqmeta.tsv | head -3 | tr '\n' ' ')"
    echo "  -> HAP1 points at a different assembly than ACC_HAP1, or the fasta was renamed."
    exit 1
  fi

  # class: assembled molecules are chromosomes, everything else is unplaced
  awk -F'\t' -v OFS="\t" 'BEGIN{print "seqid","chr_name","role","length","seq_name","class"}
      {cl = ($3 ~ /assembled/) ? "auto" : "unplaced"; print $1,$2,$3,$4,$5,cl}' tables/hap1_seqmeta.tsv > tables/chrom_classes.tsv

  # resolve a config value that may be a sequence ID, a chr_name, a seq_name (SUPER_n), or auto
  resolve_chrom() {   # $1 value  $2 label
    local v=$1 lab=$2 hit
    [[ -z $v ]] && { echo ""; return; }
    if [[ $v == auto ]]; then
      hit=$(awk -F'\t' 'NR>1 && $6=="auto"{print $4"\t"$1}' tables/chrom_classes.tsv | sort -k1,1n | head -1 | cut -f2)
    elif awk -F'\t' -v v="$v" 'NR>1 && $1==v{f=1} END{exit !f}' tables/chrom_classes.tsv; then
      hit=$v
    else
      hit=$(awk -F'\t' -v v="$v" 'NR>1 && ($2==v || $5==v){print $1}' tables/chrom_classes.tsv | head -1)
    fi
    [[ -n $hit ]] || { echo "ERROR: could not resolve $lab='$v' to a sequence in tables/chrom_classes.tsv" > /dev/stderr; exit 1; }
    echo "$hit"
  }
  MICRO=$(resolve_chrom "$MICRO" MICRO)
  NEOX=$(resolve_chrom "$NEOX" NEOX)
  NEOY=$(resolve_chrom "${NEOY:-}" NEOY)
  [[ -n $NEOX ]] || echo "WARNING: NEOX unset - the X will sit in the autosomal baseline"
  [[ -n $NEOY ]] || echo "WARNING: NEOY unset - if this assembly has a male-limited chromosome it will sit in the autosomal baseline and bias every enrichment downward"
  EXCL_IDS=""
  for e in ${EXCLUDE:-}; do EXCL_IDS="$EXCL_IDS $(resolve_chrom "$e" EXCLUDE)"; done

  # apply the special classes (exclusions last so they always win)
  awk -F'\t' -v OFS="\t" -v mi="$MICRO" -v nx="$NEOX" -v ny="$NEOY" -v ex="$EXCL_IDS" \
      'BEGIN{n=split(ex,E," "); for(i=1;i<=n;i++) X[E[i]]=1}
       NR==1{print; next}
       {if($1==mi) $6="micro"; else if(nx!="" && $1==nx) $6="neoX"; else if(ny!="" && $1==ny) $6="neoY";
        if($1 in X) $6="excluded"; print}' \
      tables/chrom_classes.tsv > tables/tmp && mv tables/tmp tables/chrom_classes.tsv
  printf 'MICRO=%s\nNEOX=%s\nNEOY=%s\nNEOX_BREAK=%s\nNEOX_ARM=%s\n' "$MICRO" "$NEOX" "$NEOY" "${NEOX_BREAK:-}" "$NEOX_ARM" > tables/resolved.env

  # print the classes and the X arm lengths
  echo "--- assembled molecules (seqid chr_name seq_name length class):"
  awk -F'\t' -v OFS="\t" 'NR==1 || $6!="unplaced"{print $1,$2,$5,$4,$6}' tables/chrom_classes.tsv | { column -t 2>/dev/null || cat; }
  echo "--- class counts:"; awk -F'\t' 'NR>1{c[$6]++} END{for(k in c) print k, c[k]}' tables/chrom_classes.tsv
  if [[ -n $NEOX && -n ${NEOX_BREAK:-} ]]; then
    xlen=$(awk -F'\t' -v x="$NEOX" '$1==x{print $4}' tables/chrom_classes.tsv)
    awk -v b="$NEOX_BREAK" -v L="$xlen" -v arm="$NEOX_ARM" 'BEGIN{
      if(b+0<=0 || b+0>=L+0){print "ERROR: NEOX_BREAK "b" is outside the NEOX scaffold (length "L")" > "/dev/stderr"; exit 1}
      la=(arm=="left")?b:L-b; ra=L-la;
      printf "neo-X arm %.1f Mb / ancestral X arm %.1f Mb (breakpoint %.3f Mb, neo-X on the %s)\n", la/1e6, ra/1e6, b/1e6, arm}'
  fi

  # female assembly metadata, for reference only (not chromosome-level)
  $DATASETS summary genome accession $ACC_FEM --report sequence --as-json-lines > tables/fem_seqreport.jsonl || true
fi
# later steps rerun standalone pick up the resolved IDs
[[ -f tables/resolved.env ]] && source tables/resolved.env
echo "resolved: MICRO=${MICRO:-<unset>} NEOX=${NEOX:-<unset>} NEOY=${NEOY:-<unset>}"

# step 1: count k-mers in reads and assemblies; solid-k-mer thresholds from the read histograms
if run_step 1; then
  meryl count k=$K memory=$MEM_GB threads=$THREADS $MALE_READS output meryl/male_reads.meryl
  meryl count k=$K memory=$MEM_GB threads=$THREADS $FEM_READS  output meryl/fem_reads.meryl
  meryl count k=$K memory=$MEM_GB threads=$THREADS $HAP1 output meryl/hap1.meryl
  meryl count k=$K memory=$MEM_GB threads=$THREADS $HAP2 output meryl/hap2.meryl
  meryl count k=$K memory=$MEM_GB threads=$THREADS $FEM  output meryl/fem_asm.meryl
  meryl histogram meryl/male_reads.meryl > tables/male_reads.hist
  meryl histogram meryl/fem_reads.meryl  > tables/fem_reads.hist
  printf "sex\tvalley\tpeak\tmin_count\n" > tables/thresholds.tsv
  for s in male fem; do
    read -r v p < <(hist_valley_peak tables/${s}_reads.hist)
    [[ "$v" == 0 ]] && { echo "no valley found in ${s} histogram - inspect tables/${s}_reads.hist and set MIN_${s^^} manually"; exit 1; }
    m=$v
    [[ $s == male && $MIN_MALE != auto ]] && m=$MIN_MALE
    [[ $s == fem  && $MIN_FEM  != auto ]] && m=$MIN_FEM
    printf "%s\t%s\t%s\t%s\n" $s $v $p $m >> tables/thresholds.tsv
    # hemizygous k-mers sit near peak/2 in the male: min_count must be well below that
    echo "$s: valley=$v peak=$p min_count=$m  (half-peak = $((p/2)))"
  done
fi

# step 3: set operations, set sizes and same-individual dropout
if run_step 3; then
  # solid k-mers per sex
  MM=$(awk '$1=="male"{print $4}' tables/thresholds.tsv); MF=$(awk '$1=="fem"{print $4}' tables/thresholds.tsv)
  meryl greater-than $((MM-1)) meryl/male_reads.meryl output meryl/male_solid.meryl
  meryl greater-than $((MF-1)) meryl/fem_reads.meryl  output meryl/fem_solid.meryl
  # strict on presence, lenient on absence: subtract the FULL opposite-sex read set, singletons included
  meryl difference meryl/male_solid.meryl meryl/fem_reads.meryl  output meryl/male_specific.meryl
  meryl difference meryl/fem_solid.meryl  meryl/male_reads.meryl output meryl/fem_specific.meryl
  # single-copy is defined from the assembly, not from read counts (hemizygous k-mers sit at half coverage in the male)
  meryl equal-to 1 meryl/hap1.meryl    output meryl/hap1_1cp.meryl
  meryl equal-to 1 meryl/fem_asm.meryl output meryl/fem_1cp.meryl
  # male tracks on hap1: denominator = read-supported hap1 k-mers; numerator = those never seen in the female
  meryl intersect meryl/hap1_1cp.meryl meryl/male_solid.meryl    output meryl/den1.meryl
  meryl intersect meryl/den1.meryl     meryl/male_specific.meryl output meryl/num1.meryl
  # repeat-inclusive track: all hap1 k-mers, not only single-copy ones
  meryl intersect meryl/hap1.meryl meryl/male_solid.meryl    output meryl/denA.meryl
  meryl intersect meryl/denA.meryl meryl/male_specific.meryl output meryl/numA.meryl
  # reciprocal control on female contigs
  meryl intersect meryl/fem_1cp.meryl meryl/fem_solid.meryl    output meryl/denF.meryl
  meryl intersect meryl/denF.meryl    meryl/fem_specific.meryl output meryl/numF.meryl
  # where male-specific k-mers sit: hap1 / hap2 only / in neither assembly
  meryl intersect  meryl/male_specific.meryl meryl/hap1.meryl output meryl/ms_hap1.meryl
  meryl difference meryl/male_specific.meryl meryl/hap1.meryl output meryl/ms_nothap1.meryl
  meryl intersect  meryl/ms_nothap1.meryl    meryl/hap2.meryl output meryl/ms_hap2only.meryl
  # same-individual dropout: single-copy assembly k-mers never seen in that individual's own reads at any count
  meryl difference meryl/hap1_1cp.meryl meryl/male_reads.meryl output meryl/x_hap1_1cp_unseen.meryl
  meryl difference meryl/fem_1cp.meryl  meryl/fem_reads.meryl  output meryl/x_fem_1cp_unseen.meryl
  # set sizes and dropout tables
  printf "set\tdistinct_kmers\n" > tables/set_sizes.tsv
  for s in male_reads fem_reads male_solid fem_solid male_specific fem_specific hap1 hap2 fem_asm \
           hap1_1cp fem_1cp den1 num1 denF numF ms_hap1 ms_hap2only ms_nothap1; do
    printf "%s\t%s\n" $s "$(nk meryl/$s.meryl)" >> tables/set_sizes.tsv
  done
  a=$(nk meryl/hap1_1cp.meryl); b=$(nk meryl/x_hap1_1cp_unseen.meryl)
  c=$(nk meryl/fem_1cp.meryl);  d=$(nk meryl/x_fem_1cp_unseen.meryl)
  awk -v a=$a -v b=$b -v c=$c -v d=$d -v OFS="\t" \
    'BEGIN{print "individual","asm_1cp_kmers","unseen_in_own_reads","dropout"; print "male",a,b,b/a; print "fem",c,d,d/c}' > tables/dropout.tsv
  cat tables/set_sizes.tsv tables/dropout.tsv
fi

# step 4: place k-mers on hap1 (male tracks) and on female contigs (reciprocal); window counts
if run_step 4; then
  [[ -s tables/hap1.genome ]] || seqkit fx2tab -n -i -l $HAP1 > tables/hap1.genome
  seqkit fx2tab -n -i -l $FEM  > tables/fem.genome
  lookup_windows $HAP1 meryl/den1.meryl den1 meryl/num1.meryl num1
  lookup_windows $HAP1 meryl/denA.meryl denA meryl/numA.meryl numA
  lookup_windows $FEM meryl/denF.meryl denF meryl/numF.meryl numF
  # browsable BED of single-copy male-specific k-mers on hap1 (adjacent k-mers merged)
  meryl-lookup -bed-runs -sequence $HAP1 -mers meryl/num1.meryl -output beds/male_specific_1cp_on_hap1.bed
fi

# step 5: read depth of both sexes on hap1 (mosdepth, W-bp windows, MAPQ>=20)
if run_step 5; then
  for s in male fem; do
    if [[ $s == male ]]; then R="$MALE_READS"; P=$MM_PRESET_MALE; else R="$FEM_READS"; P=$MM_PRESET_FEM; fi
    minimap2 -t $THREADS -ax $P --secondary=no $HAP1 $R \
    | samtools sort -@ 8 -m 3G -T aln/tmp_$s -o aln/${s}_hap1.bam -
    samtools index aln/${s}_hap1.bam
    mosdepth -t 8 -n -Q 20 --by $W aln/${s}_hap1 aln/${s}_hap1.bam
  done
fi

# step 6: assembly alignments and coverage tables
if run_step 6; then
  # female contigs -> hap1 (asm20 tolerates cross-individual divergence; secondaries kept to see diverged homology)
  minimap2 -t $THREADS -cx asm20 --secondary=yes -N 20 $HAP1 $FEM > aln/fem_to_hap1.paf
  # hap2 -> hap1 at asm5: only near-identical (allelic) sequence aligns; hemizygous chromosomes should be ~uncovered
  minimap2 -t $THREADS -cx asm5 --secondary=no $HAP1 $HAP2 > aln/hap2_to_hap1.paf
  paf_cov aln/fem_to_hap1.paf   P   target tables/hap1.genome  tables/asmcov_fem_on_hap1_primary.tsv
  paf_cov aln/fem_to_hap1.paf   ALL target tables/hap1.genome  tables/asmcov_fem_on_hap1_all.tsv
  paf_cov aln/hap2_to_hap1.paf  P   target tables/hap1.genome  tables/asmcov_hap2_on_hap1_primary.tsv
  paf_slim aln/fem_to_hap1.paf aln/fem_to_hap1.slim.tsv
  # each candidate sex-limited chromosome vs the rest of hap1 (diverged homolog) and vs hap2 (allelic partner)
  for role in micro neoy; do
    [[ $role == micro ]] && cid=$MICRO || cid=${NEOY:-}
    [[ -n $cid ]] || continue
    seqkit grep -p  $cid $HAP1 > aln/$role.fa
    seqkit grep -v -p $cid $HAP1 > aln/hap1_no_$role.fa
    seqkit fx2tab -n -i -l aln/$role.fa > tables/$role.genome
    minimap2 -t $THREADS -cx asm20 --secondary=yes -N 20 aln/hap1_no_$role.fa aln/$role.fa > aln/${role}_to_hap1rest.paf
    minimap2 -t $THREADS -cx asm20 --secondary=yes -N 20 $HAP2 aln/$role.fa > aln/${role}_to_hap2.paf
    paf_cov aln/${role}_to_hap1rest.paf ALL query tables/$role.genome tables/asmcov_${role}_by_hap1rest_all.tsv
    paf_cov aln/${role}_to_hap2.paf     ALL query tables/$role.genome tables/asmcov_${role}_by_hap2_all.tsv
    paf_slim aln/${role}_to_hap1rest.paf aln/${role}_to_hap1rest.slim.tsv
    echo "$role ($cid) covered by rest of hap1 / by hap2 (all alignments):"
    cat tables/asmcov_${role}_by_hap1rest_all.tsv tables/asmcov_${role}_by_hap2_all.tsv
    rm -f aln/hap1_no_$role.fa
  done
fi

# step 7: per-sequence length / GC (light contamination screen)
if run_step 7; then
  seqkit fx2tab -n -i -l -g $HAP1 > tables/hap1_seqstats.tsv
fi

# step 8: R summary
if run_step 8; then
  Rscript "$SCRIPT_DIR/summarize_kmer_sex.R" . $W "$MICRO" "$NEOX" $MIN_CHR_LEN $NPERM "${NEOX_BREAK:-NA}" "$NEOX_ARM" "${NEOY:-}" | tee tables/summary.txt
fi

# step 9: heterozygosity per class from the step-5 BAMs (hemizygous vs mosaic copy number)
# slow: best submitted alone with STEPS="9"; it re-runs the R summary at the end
if run_step 9; then
  command -v bcftools >/dev/null || { echo "ERROR: bcftools not found - conda install -n kmer -c bioconda bcftools"; exit 1; }
  bcftools --version | head -1
  ln -sf "$HAP1" aln/hap1.fa; [[ -s aln/hap1.fa.fai ]] || samtools faidx aln/hap1.fa
  mkdir -p aln/hetlog
  # regions: the candidates plus the N_HET_AUTO largest autosomes (keeps mpileup to a few hundred Mb)
  REGS="$(awk -F'\t' 'NR>1 && $6=="auto"{print $4"\t"$1}' tables/chrom_classes.tsv | sort -k1,1nr | head -$N_HET_AUTO | cut -f2 | tr '\n' ' ')"
  REGS="$REGS $MICRO ${NEOX:-} ${NEOY:-}"
  printf '%s\n' $REGS > tables/het_regions.txt
  # one region is enough to prove the command works; fail here rather than 6 hours in
  probe=$(head -1 tables/het_regions.txt)
  bcftools mpileup -f aln/hap1.fa -r "$probe:1-100000" -q 20 -Q 20 -d 500 -I \
      -a FORMAT/AD,FORMAT/DP aln/male_hap1.bam 2> aln/hetlog/probe.log \
    | bcftools call -mv -Ov > /dev/null 2>> aln/hetlog/probe.log \
    || { echo "ERROR: bcftools probe on $probe failed:"; cat aln/hetlog/probe.log; exit 1; }
  echo "bcftools probe OK on $probe"

  NPAR=${NPAR:-6}                       # each worker streams the whole BAM region; more is not faster on Lustre
  for s in male fem; do
    # call SNPs per region in parallel (bash -c + pipefail so a failed mpileup is not masked)
    printf '%s\n' $REGS | xargs -P $NPAR -I{} bash -c \
      "set -o pipefail; bcftools mpileup -f aln/hap1.fa -r {} -q 20 -Q 20 -d 500 -I -a FORMAT/AD,FORMAT/DP aln/${s}_hap1.bam 2> aln/hetlog/${s}_{}.log \
       | bcftools call -mv -Oz -o aln/${s}_het_{}.vcf.gz 2>> aln/hetlog/${s}_{}.log && bcftools index -f aln/${s}_het_{}.vcf.gz" \
      || { echo "ERROR: step 9 calling failed for $s. Region logs:"; tail -n 5 aln/hetlog/${s}_*.log; exit 1; }
    bcftools concat -a -Oz -o aln/${s}_hap1.calls.vcf.gz $(for r in $REGS; do echo aln/${s}_het_$r.vcf.gz; done)
    bcftools index -f aln/${s}_hap1.calls.vcf.gz
    # numerator: het SNPs inside the depth band -> per window n_het, mean alt fraction, n balanced (0.3-0.7)
    bcftools view -v snps -g het -i "QUAL>=30 && FMT/DP>=$HET_DP_MIN && FMT/DP<=$HET_DP_MAX" aln/${s}_hap1.calls.vcf.gz \
      | bcftools query -f '%CHROM\t%POS\t[%AD]\n' \
      | awk -F'\t' -v W=$W -v OFS="\t" '{split($3,a,","); if(a[1]+a[2]==0) next; af=a[2]/(a[1]+a[2]);
            b=int(($2-1)/W); k=$1 SUBSEP b; n[k]++; sum[k]+=af; if(af>=0.3 && af<=0.7) bal[k]++}
         END{for(k in n){split(k,x,SUBSEP); print x[1], x[2]*W, (x[2]+1)*W, n[k], sum[k]/n[k], bal[k]+0}}' \
      | sort -k1,1 -k2,2n > tables/${s}_het_windows.tsv
    [[ -s tables/${s}_het_windows.tsv ]] || echo "WARNING: no het SNPs passed filters for $s"
    # denominator: callable sites per window in the same depth band and base/mapping quality
    printf '%s\n' $REGS | xargs -P $NPAR -I{} bash -c \
      "set -o pipefail; samtools depth -q 20 -Q 20 -r {} aln/${s}_hap1.bam \
       | awk -v W=$W -v lo=$HET_DP_MIN -v hi=$HET_DP_MAX -v OFS='\t' \
           '\$3>=lo && \$3<=hi {b=int((\$2-1)/W); c[\$1 SUBSEP b]++}
            END{for(k in c){split(k,x,SUBSEP); print x[1], x[2]*W, (x[2]+1)*W, c[k]}}' \
       > aln/hetlog/${s}_callable_{}.tsv" \
      || { echo "ERROR: step 9 callable-site pass failed for $s"; exit 1; }
    cat aln/hetlog/${s}_callable_*.tsv | sort -k1,1 -k2,2n > tables/${s}_callable_windows.tsv
    rm -f aln/hetlog/${s}_callable_*.tsv
    echo "$s callable sites in DP band [$HET_DP_MIN,$HET_DP_MAX] per region:"
    awk -F'\t' '{c[$1]+=$4} END{for(k in c) printf "  %s\t%d\n", k, c[k]}' tables/${s}_callable_windows.tsv
    rm -f aln/${s}_het_*.vcf.gz aln/${s}_het_*.vcf.gz.csi
  done
  # het results feed the summary, so refresh it
  Rscript "$SCRIPT_DIR/summarize_kmer_sex.R" . $W "$MICRO" "$NEOX" $MIN_CHR_LEN $NPERM "${NEOX_BREAK:-NA}" "$NEOX_ARM" "${NEOY:-}" | tee tables/summary_with_het.txt
fi
echo done
