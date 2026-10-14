#!/bin/bash
#SBATCH -D /path/to/project/gene_pairs_ds
#SBATCH -J ds_multi
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH -o /path/to/project/gene_pairs_ds/slurm-%j.out
#SBATCH -e /path/to/project/gene_pairs_ds/slurm-%j.err
# Gene pairs, dN/dS and the fusion breakpoint for A. juncea fused X / neo-Y vs Anax junius X / chr12 (Fig. 3C)

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK:-16}"
WORKDIR="$(pwd)"
C="$WORKDIR/chrom"
GUIDANCE_PY="${GUIDANCE_PY:-/path/to/project/programs/guidance_pairwise.py}"
DATA_GENOMES="/path/to/project/data/genomes"      # <species>$GENOME_SUFFIX
DATA_GENES="/path/to/project/data/genes"          # <species>$GFF_SUFFIX

# suffixes of the data-drop inputs (the chrom/ files written here keep plain .fna/.gff3 names);
# the genomes must be soft-masked: hard-masked N inside CDS would corrupt every dN/dS
GENOME_SUFFIX=".clean.masked.fna"
GFF_SUFFIX=".longest.gff"

MINIPROT="${MINIPROT:-miniprot}"       # set full path if not on PATH

# entity label -> species (locates $DATA_GENOMES/<sp>$GENOME_SUFFIX + $DATA_GENES/<sp>$GFF_SUFFIX) + scaffold
declare -A SPECIES SCAFFOLD
SPECIES[juncea_X]="Aeshna.juncea";       SCAFFOLD[juncea_X]="OZ369554.1"   # fused X
SPECIES[juncea_Y]="Aeshna.juncea";       SCAFFOLD[juncea_Y]="OZ369555.1"   # neo-Y
SPECIES[anax_X]="Anax.junius";           SCAFFOLD[anax_X]="OZ349813.1"     # Anax (ancestral) X
SPECIES[anax_12]="Anax.junius";          SCAFFOLD[anax_12]="OZ349811.1"    # Anax chr12 = unfused neo-X homolog
SPECIES[soma_X]="Somatochlora.arctica";  SCAFFOLD[soma_X]="OZ407525.1"     # Somatochlora ancestral X (built in the v8 run; not used downstream)

# comparisons "A B"; the pair files are oriented geneA<TAB>geneB
COMPARISONS=(
  "juncea_X juncea_Y"     # gametologs
  "juncea_X anax_12"      # neo-X vs Anax 12
  "juncea_Y anax_12"      # neo-Y vs Anax 12
  "juncea_X anax_X"       # ancestral X vs Anax X
)

# activate the env with set -u/-x off: the deactivate hook of the env active at submission reads an
# unset variable (fatal under nounset, before logging starts), and -x would print hundreds of lines
{ set +ux; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/ds   # agat,gffread,samtools,diamond,MCScanX,orthofinder,codeml,miniprot
set -ux

# preflight: the env and binaries are checked now rather than hours in
[[ "${CONDA_PREFIX:-}" == "/path/to/project/conda/envs/ds" ]] \
  || { echo "[error] conda env not active: CONDA_PREFIX=${CONDA_PREFIX:-unset}" >&2; exit 1; }

NEED_BIN=(python3 diamond MCScanX orthofinder codeml gffread samtools
          agat_sp_keep_longest_isoform.pl "$MINIPROT")
for b in "${NEED_BIN[@]}"; do
  command -v "$b" >/dev/null || { echo "[error] $b not on PATH after activate" >&2; exit 1; }
done

# biopython (process_pair converts the filtered alignment with Bio.SeqIO) and the guidance script
python3 -c 'import Bio' 2>/dev/null \
  || { echo "[error] biopython missing: conda install -n ds -c conda-forge biopython" >&2; exit 1; }
[[ -s "$GUIDANCE_PY" ]] \
  || { echo "[error] guidance script not found: $GUIDANCE_PY" >&2; exit 1; }
# OrthoFinder 2.5.5 fails on numpy >= 2 / scipy >= 1.11 and misreports it as running out of RAM
python3 - <<PYCHK || exit 1
import sys
try:
    import numpy, scipy
except Exception as e:
    sys.exit("[error] numpy/scipy missing: %s" % e)
nv = tuple(int(x) for x in numpy.__version__.split(".")[:2])
sv = tuple(int(x) for x in scipy.__version__.split(".")[:2])
if nv >= (2, 0) or sv >= (1, 11):
    sys.exit("[error] numpy %s / scipy %s: OrthoFinder 2.5.5 will fail here and "
             "report it as running out of RAM. Pin numpy<2 and scipy<1.11."
             % (numpy.__version__, scipy.__version__))
print("[preflight] numpy %s, scipy %s, python %d.%d -- OK for OrthoFinder 2.5.5"
      % (numpy.__version__, scipy.__version__,
         sys.version_info[0], sys.version_info[1]))
PYCHK

# MCScanX's recursive chaining overflows the default stack on the neo-X/neo-Y block: raise it to a large
# finite value, not `unlimited` (that switches to the legacy mmap layout, where MCScanX segfaults)
ulimit -s 524288 2>/dev/null || ulimit -s unlimited 2>/dev/null || true
echo "[preflight] stack limit: $(ulimit -s) KB"

# pin MCScanX to the active env's binary so no PATH change can shadow it
MCSCANX_BIN="${MCSCANX_BIN:-$CONDA_PREFIX/bin/MCScanX}"
[[ -x "$MCSCANX_BIN" ]] \
  || { echo "[error] MCScanX not found at $MCSCANX_BIN" >&2; exit 1; }
echo "[preflight] MCScanX: $MCSCANX_BIN"
echo "[preflight] first MCScanX on PATH: $(command -v MCScanX || echo none)"

# output folders, log to logs/, temporary files in tmp/
mkdir -p chrom comp fate logs tmp
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)
TMPROOT="$WORKDIR/tmp"
export TMPDIR="$TMPROOT"

# an entity's source GFF and genome
gff_for()    { echo "$DATA_GENES/${SPECIES[$1]}${GFF_SUFFIX}"; }
genome_for() { echo "$DATA_GENOMES/${SPECIES[$1]}${GENOME_SUFFIX}"; }

# per-chromosome build (cached in chrom/): scaffold FASTA, longest-isoform GFF, CDS, protein, ids
build_entity() {
  local L="$1"
  local FNA="$C/${L}.fna"      GFF3="$C/${L}.gff3"  LGFF="$C/${L}.longest.gff3"
  local CDS="$C/${L}.cds.fna"  PEP="$C/${L}.pep.faa" IDS="$C/${L}.ids"

  if [[ -s "$CDS" && -s "$CDS.fai" && -s "$PEP" && -s "$IDS" && -s "$LGFF" ]]; then
    echo "[chrom] $L already built, skipping"; return 0
  fi

  local SCAF="${SCAFFOLD[$L]}" SRCG SRCGFF
  SRCG="$(genome_for "$L")"; SRCGFF="$(gff_for "$L")"
  [[ -f "$SRCG"   ]] || { echo "[error] genome not found: $SRCG" >&2; exit 1; }
  [[ -f "$SRCGFF" ]] || { echo "[error] gff not found: $SRCGFF" >&2; exit 1; }

  echo "[chrom] building $L  (${SPECIES[$L]} : $SCAF)"
  samtools faidx "$SRCG" "$SCAF" > "$FNA"
  [[ -s "$FNA" ]] || { echo "[error] $SCAF not found in $SRCG" >&2; exit 1; }

  awk -F'\t' -v c="$SCAF" '$1==c' "$SRCGFF" > "$GFF3"
  [[ -s "$GFF3" ]] || { echo "[error] no features for $SCAF in $SRCGFF" >&2; exit 1; }

  # longest isoform per gene (a no-op on the .longest.gff inputs; kept so the step does not depend on them)
  agat_sp_keep_longest_isoform.pl --gff "$GFF3" -o "$LGFF"
  gffread "$LGFF" -g "$FNA" -x "$CDS" -y "$PEP"
  grep '^>' "$PEP" | sed 's/^>//; s/[[:space:]].*//' | sort -u > "$IDS"
  samtools faidx "$CDS"
}

# uncapped (-k 0) DIAMOND search for the anchor-partner ranks in lib/anchor_filter.R (ab.blast stays -k 5
# for MCScanX); a separate file, not filtered, built even when the union is cached
build_rank_blast() {
  local A="$1" B="$2" CMP="$3"
  local C="$WORKDIR/chrom"
  [[ -s "$CMP/ab.rank.blast" ]] && { echo "[rank] $CMP exists, skipping"; return 0; }

  echo "[rank] deep search for $A vs $B (-k 0, uncapped)"
  [[ -s "$CMP/ab.pep.faa" ]] || cat "$C/${A}.pep.faa" "$C/${B}.pep.faa" > "$CMP/ab.pep.faa"
  [[ -s "$CMP/ab.dmnd" ]]    || diamond makedb --in "$CMP/ab.pep.faa" -d "$CMP/ab.dmnd"

  diamond blastp -q "$CMP/ab.pep.faa" -d "$CMP/ab.dmnd" -o "$CMP/ab.rank.blast.tmp" \
    --outfmt 6 -e 1e-5 -k 0 -p "$THREADS" --more-sensitive
  sort -k1,1 -k2,2 "$CMP/ab.rank.blast.tmp" > "$CMP/ab.rank.blast"
  rm -f "$CMP/ab.rank.blast.tmp"

  # hit-count summary (unique query-subject pairs, self-hit included); a pile-up at the max means a cap
  awk -F'\t' '!seen[$1 FS $2]++ {n[$1]++} END {
    mx=0; for (q in n) if (n[q] > mx) mx = n[q]
    c=0;  for (q in n) if (n[q] == mx) c++
    printf "[rank] %s: %d queries, max %d hits, %d queries at that max\n",
           "'"$CMP"'", length(n), mx, c
  }' "$CMP/ab.rank.blast"
}

# per-comparison ortholog union (cached as comp/<name>/pairs.union.tsv): DIAMOND + MCScanX + OrthoFinder
build_union() {
  local A="$1" B="$2" CMP="$3"
  local C="$WORKDIR/chrom"
  local PAIRS="$CMP/pairs.union.tsv"
  [[ -s "$PAIRS" ]] && { echo "[union] $CMP exists, skipping"; return 0; }

  echo "[union] $A vs $B"#

  # diamond all-vs-all
  cat "$C/${A}.pep.faa" "$C/${B}.pep.faa" > "$CMP/ab.pep.faa"
  diamond makedb --in "$CMP/ab.pep.faa" -d "$CMP/ab.dmnd"
  diamond blastp -q "$CMP/ab.pep.faa" -d "$CMP/ab.dmnd" -o "$CMP/ab.blast" \
    --outfmt 6 -e 1e-5 -k 5 -p "$THREADS" --more-sensitive

  # MCScanX gff: the two chromosomes tagged as "genomes" aa and bb
  awk -v T=aa '$3=="mRNA"{match($9,/ID=[^;]+/); id=substr($9,RSTART+3,RLENGTH-3); print T"\t"id"\t"$4"\t"$5}' \
      "$C/${A}.longest.gff3" >  "$CMP/ab.gff"
  awk -v T=bb '$3=="mRNA"{match($9,/ID=[^;]+/); id=substr($9,RSTART+3,RLENGTH-3); print T"\t"id"\t"$4"\t"$5}' \
      "$C/${B}.longest.gff3" >> "$CMP/ab.gff"

  # keep only hits with both genes in ab.gff (any other makes MCScanX segfault), sorted for determinism
  cut -f2 "$CMP/ab.gff" | sort -u > "$CMP/ab.genes"
  awk 'NR==FNR{g[$1]; next} ($1 in g) && ($2 in g)' "$CMP/ab.genes" "$CMP/ab.blast" \
      | sort -k1,1 -k2,2 > "$CMP/ab.blast.tmp" && mv "$CMP/ab.blast.tmp" "$CMP/ab.blast"

  # run MCScanX, output to mcscanx.log; a non-zero exit is only a warning, since its HTML step
  # (after the collinearity is written) can segfault in batch
  if ! ( cd "$CMP" && "$MCSCANX_BIN" -b 2 ab > mcscanx.log 2>&1 ); then
    echo "[warn] MCScanX exited non-zero for $CMP -- its output:" >&2
    sed 's/^/[mcscanx] /' "$CMP/mcscanx.log" >&2 || true
  fi

  # its anchor pairs -> mcscanx.pairs.tsv; no collinearity is fatal (figs/fig3/fig3c.R would get no anchors)
  [[ -s "$CMP/ab.collinearity" ]] || {
    echo "[error] MCScanX produced no collinearity for $CMP." >&2
    echo "[error] It crashed during chaining, not in the HTML step. Reproduce with:" >&2
    echo "[error]   cd $CMP && $MCSCANX_BIN -b 2 ab" >&2
    exit 1
  }
  awk '!/^#/ && NF>=5 {print $3"\t"$4}' "$CMP/ab.collinearity" | sort -u > "$CMP/mcscanx.pairs.tsv"
  [[ -s "$CMP/mcscanx.pairs.tsv" ]] || {
    echo "[error] ab.collinearity exists but yielded no pairs for $CMP -- check its format" >&2
    exit 1
  }
  echo "[union] $(wc -l < "$CMP/mcscanx.pairs.tsv") MCScanX anchor pairs"

  # orthofinder (a = side A, b = side B; fixed names keep the output path stable)
  rm -rf "$CMP/of_input" "$CMP/of_output"
  mkdir -p "$CMP/of_input"
  cp "$C/${A}.pep.faa" "$CMP/of_input/a.fa"
  cp "$C/${B}.pep.faa" "$CMP/of_input/b.fa"
  orthofinder -f "$CMP/of_input" -t "$THREADS" -o "$CMP/of_output" -S diamond_ultra_sens

  # every a-b pair of the orthologue table, and of each orthogroup
  local OF_DIR OF_TSV OG_TSV
  OF_DIR=$(ls -dt "$CMP"/of_output/Results_*/ | head -1)   # newest, not lexical-first
  OF_TSV="${OF_DIR}/Orthologues/Orthologues_a/a__v__b.tsv"
  awk -F'\t' 'NR>1 {
    n = split($2, xs, ", "); m = split($3, ys, ", ")
    for (i=1; i<=n; i++) for (j=1; j<=m; j++) print xs[i]"\t"ys[j]
  }' "$OF_TSV" | sort -u > "$CMP/orthofinder.pairs.tsv"

  OG_TSV="${OF_DIR}/Orthogroups/Orthogroups.tsv"
  awk -F'\t' 'NR>1 && $2!="" && $3!="" {
    n = split($2, xs, ", "); m = split($3, ys, ", ")
    for (i=1; i<=n; i++) for (j=1; j<=m; j++) print xs[i]"\t"ys[j]
  }' "$OG_TSV" | sort -u > "$CMP/orthogroups.pairs.tsv"

  # union of all sources, oriented as A<TAB>B, intra-chromosome pairs dropped
  cat "$CMP/mcscanx.pairs.tsv" "$CMP/orthofinder.pairs.tsv" "$CMP/orthogroups.pairs.tsv" \
  | awk -F'\t' -v af="$C/${A}.ids" -v bf="$C/${B}.ids" '
      BEGIN{ while((getline k < af)>0) AA[k]; while((getline k < bf)>0) BB[k] }
      ($1 in AA && $2 in BB){ print $1"\t"$2; next }
      ($1 in BB && $2 in AA){ print $2"\t"$1 }
    ' | sort -u > "$PAIRS"
  sed -i 's/\r//g' "$PAIRS"

  # the union must contain every MCScanX anchor (a smaller one is a stale cached pairs.union.tsv)
  local NU NM
  NU=$(wc -l < "$PAIRS"); NM=$(wc -l < "$CMP/mcscanx.pairs.tsv")
  echo "[union] $CMP: $NU union pairs (>= $NM anchors)"
  (( NU >= NM )) || {
    echo "[error] union ($NU) smaller than the MCScanX anchor set ($NM) for $CMP." >&2
    echo "[error] Delete $PAIRS and rebuild." >&2
    exit 1
  }

  # advisory per-pair coverage (not a filter, used by figs/fig3/fig3c.R): best single-HSP coverage of
  # gene_a and gene_b, and their minimum; 0 for a pair with no direct DIAMOND hit
  [[ -s "$CMP/ab.pep.faa.fai" ]] || samtools faidx "$CMP/ab.pep.faa"
  awk -F'\t' -v plen="$CMP/ab.pep.faa.fai" -v pf="$PAIRS" '
    BEGIN{
      while((getline l < plen)>0){ split(l,t,"\t"); L[t[1]]=t[2] }
      while((getline l < pf)>0){ split(l,p,"\t"); want[p[1] SUBSEP p[2]]=1 }
    }
    # diamond outfmt6 default: 1=q 2=s 7=qstart 8=qend 9=sstart 10=send
    { q=$1; s=$2
      if(L[q]>0){ c=($8-$7+1)/L[q]; if(c>cov[q SUBSEP s]) cov[q SUBSEP s]=c }
      if(L[s]>0){ c=($10-$9+1)/L[s]; if(c>cov[s SUBSEP q]) cov[s SUBSEP q]=c } }
    END{
      for(k in want){ split(k,p,SUBSEP); a=p[1]; b=p[2]
        ca=cov[a SUBSEP b]+0; cb=cov[b SUBSEP a]+0; mc=(ca<cb?ca:cb)
        printf "%s\t%s\t%.3f\t%.3f\t%.3f\n", a,b,ca,cb,mc }
    }' "$CMP/ab.blast" | sort -k1,1 -k2,2 > "$CMP/pairs.cover.tsv"
}

# per-pair GUIDANCE alignment + codeml pairwise ML; reads CMP_DIR/PAIRS_FILE/A_CDS/B_CDS from the loop
process_pair() {
  local TASK_ID="$1"
  local OUTDIR="$CMP_DIR/pairs/pair_${TASK_ID}"
  local PTMP="$TMPROOT/${CMP_NAME}_task_${TASK_ID}"
  mkdir -p "$OUTDIR" "$PTMP"

  # the pair on line TASK_ID of the current union (pair_N is keyed by line number)
  local LINE AID BID
  LINE=$(sed -n "${TASK_ID}p" "$PAIRS_FILE")
  AID=$(echo "$LINE" | cut -f1)   # side-A gene
  BID=$(echo "$LINE" | cut -f2)   # side-B gene

  # resume: skip only if the cached result is for this exact pair (gene ids checked) and has a real dS
  if [[ -s "$OUTDIR/result.tsv" ]]; then
    local cA cB cDS
    { IFS=$'\t' read -r cA cB _ cDS _ || true; } < "$OUTDIR/result.tsv"
    if [[ "$cA" == "$AID" && "$cB" == "$BID" && -n "${cDS:-}" && "${cDS:-}" != "NA" ]]; then
      rm -rf "$PTMP"; return 0
    fi
  fi

  (
    set -euo pipefail
    export TMPDIR="$PTMP"

    samtools faidx "$A_CDS" "$AID" >  "$OUTDIR/pair.fa"
    samtools faidx "$B_CDS" "$BID" >> "$OUTDIR/pair.fa"

    # alignment confidence: keep the last stdout line that is exactly two numerics
    GUIDEDIR="$OUTDIR/guidance"
    SCORE_LINE=$(
      python3 "$GUIDANCE_PY" \
        --seqfile "$OUTDIR/pair.fa" \
        --outdir "$GUIDEDIR" \
        --bootstraps 300 \
        --col_cutoff 0.93 \
      | grep -E '^[[:space:]]*[0-9.eE+-]+[[:space:]]+[0-9.eE+-]+[[:space:]]*$' | tail -1 || true
    )
    read -r MEAN_SCORE FRAC_RETAINED <<< "$SCORE_LINE"
    if [[ -z "${MEAN_SCORE:-}" || -z "${FRAC_RETAINED:-}" ]]; then
      echo "[warn] no guidance scores parsed for pair ${TASK_ID}" >&2
      MEAN_SCORE="NA"; FRAC_RETAINED="NA"
    fi

    # filtered alignment -> paml sequential phylip
    python3 -c "
from Bio import SeqIO
import sys
recs = list(SeqIO.parse(sys.argv[1], 'fasta'))
with open(sys.argv[2], 'w') as f:
    f.write(f' {len(recs)} {len(recs[0].seq)}\n')
    for r in recs:
        f.write(f'{r.id}\n{str(r.seq)}\n')
" "$GUIDEDIR/filtered.fa" "$OUTDIR/pair.phy"

    # codeml pairwise ML (runmode = -2)
    cat > "$OUTDIR/codeml.ctl" << EOF
seqfile = pair.phy
outfile = codeml.out
noisy = 0
verbose = 0
runmode = -2
seqtype = 1
CodonFreq = 2
model = 0
NSsites = 0
icode = 0
fix_kappa = 0
kappa = 2
fix_omega = 0
omega = 0.5
cleandata = 0
EOF
    ( cd "$OUTDIR" && codeml codeml.ctl )

    # dN, dS and omega from the last dN/dS line of codeml.out
    read -r DN DS OMEGA < <(
      awk '/dN\/dS/{
        for(i=1;i<=NF;i++){
          if($i=="dN/dS=")    w=$(i+1);
          else if($i=="dN")  dn=$(i+2);
          else if($i=="dS")  ds=$(i+2);
        }
        print dn, ds, w
      }' "$OUTDIR/codeml.out" | tail -1
    )
    # no dN/dS line from codeml (e.g. no usable codons) -> NA, so the row keeps its 7 fields
    : "${DN:=NA}"; : "${DS:=NA}"; : "${OMEGA:=NA}"

    echo -e "${AID}\t${BID}\t${DN}\t${DS}\t${OMEGA}\t${MEAN_SCORE}\t${FRAC_RETAINED}" \
      > "$OUTDIR/result.tsv"
  ) > "$OUTDIR/pair.log" 2>&1 || echo "[warn] pair ${TASK_ID} failed (see $OUTDIR/pair.log)"

  # one row per union pair: a pair that failed before writing result.tsv gets an NA row (retried next run)
  [[ -s "$OUTDIR/result.tsv" ]] \
    || printf '%s\t%s\tNA\tNA\tNA\tNA\tNA\n' "$AID" "$BID" > "$OUTDIR/result.tsv"

  rm -rf "$PTMP"; return 0
}
export -f process_pair
export GUIDANCE_PY TMPROOT

# miniprot the query proteins of a scaffold that have no pair in PRS against a target (-> fate/<name>/probe.gff)
probe_fate() {
  local NM="$1" QPEP="$2" QGFF="$3" TFNA="$4" PRS="$5" KCOL="$6" QSCAF="$7" B="$8" SIDE="$9"
  local f
  for f in "$QPEP" "$QGFF" "$TFNA"; do
    [[ -s "$f" ]] || { echo "[err] $NM: missing input $f" >&2; exit 1; }
  done
  [[ -f "$PRS" ]] || { echo "[err] $NM: pair table not found $PRS" >&2; exit 1; }   # PRS may be empty (probe all)
  local OD="fate/$NM"; mkdir -p "$OD"

  # gene universe = query mRNA ids on QSCAF, on the requested side of position B (or all)
  awk -F'\t' -v c="$QSCAF" -v b="$B" -v side="$SIDE" '
    $1==c && $3=="mRNA"{
      match($9,/ID=[^;]+/); id=substr($9,RSTART+3,RLENGTH-3)
      if      (side=="all")             print id
      else if (side=="above" && $4 > b) print id
      else if (side=="below" && $5 < b) print id
    }' "$QGFF" | sort -u > "$OD/universe.ids"

  # already-paired = query ids present in the ortholog table's query-side column
  cut -f"$KCOL" "$PRS" | sort -u > "$OD/paired.ids"
  comm -23 "$OD/universe.ids" "$OD/paired.ids" > "$OD/unpaired.ids"

  local nu np nq
  nu=$(wc -l < "$OD/universe.ids")
  np=$(comm -12 "$OD/universe.ids" "$OD/paired.ids" | wc -l)
  nq=$(wc -l < "$OD/unpaired.ids")
  echo "[info] $NM: universe=$nu  annotated-pair=$np  to-probe=$nq"
  [[ "$nq" -gt 0 ]] || { : > "$OD/probe.gff"; echo "[info] $NM: nothing to probe"; return 0; }

  [[ -s "$QPEP.fai" ]] || samtools faidx "$QPEP"
  xargs samtools faidx "$QPEP" < "$OD/unpaired.ids" > "$OD/unpaired.faa"
  "$MINIPROT" --gff -t "$THREADS" --outn 1 "$TFNA" "$OD/unpaired.faa" > "$OD/probe.gff"
  echo "[info] $NM: miniprot hits=$(grep -c '^##PAF' "$OD/probe.gff")"
}

# build every chromosome referenced by a comparison, plus soma_X
declare -A NEEDED
for cmp in "${COMPARISONS[@]}"; do read -r A B <<< "$cmp"; NEEDED[$A]=1; NEEDED[$B]=1; done
for L in "${!NEEDED[@]}"; do build_entity "$L"; done
build_entity soma_X                       # built in the v8 run; nothing downstream reads it

# headers of the combined dN/dS and coverage tables
COMBINED="$WORKDIR/all_ds_results.tsv"
echo -e "comparison\tgene_a\tgene_b\tdN\tdS\tomega\tmean_guidance_score\tfrac_cols_retained" > "$COMBINED"
COVALL="$WORKDIR/all_pairs_cover.tsv"
echo -e "comparison\tgene_a\tgene_b\tcov_a\tcov_b\tmin_cov" > "$COVALL"

# per comparison: ortholog union, rank search, dN/dS of every pair in parallel, then the tables
for cmp in "${COMPARISONS[@]}"; do
  read -r A B <<< "$cmp"
  NAME="${A}__vs__${B}"
  CMP="$WORKDIR/comp/$NAME"
  mkdir -p "$CMP/pairs"

  build_union "$A" "$B" "$CMP"
  build_rank_blast "$A" "$B" "$CMP"

  PAIRS="$CMP/pairs.union.tsv"
  NP=$(wc -l < "$PAIRS")
  echo "[info] $NAME : $NP pairs across $THREADS cores"

  # per-comparison context for process_pair
  export CMP_DIR="$CMP" CMP_NAME="$NAME" PAIRS_FILE="$PAIRS"
  export A_CDS="$WORKDIR/chrom/${A}.cds.fna" B_CDS="$WORKDIR/chrom/${B}.cds.fna"

  if (( NP > 0 )); then
    seq 1 "$NP" | xargs -P "$THREADS" -I {} bash -c 'process_pair "$@"' _ {}
  fi

  # per-comparison table (gene columns labelled by entity)
  echo -e "${A}\t${B}\tdN\tdS\tomega\tmean_guidance_score\tfrac_cols_retained" > "$CMP/ds_results.tsv"
  cat "$CMP"/pairs/pair_*/result.tsv >> "$CMP/ds_results.tsv" 2>/dev/null || true
  DONE=$(( $(wc -l < "$CMP/ds_results.tsv") - 1 ))
  echo "[done] $NAME : $DONE / $NP pairs"

  # append to the combined long-format table
  tail -n +2 "$CMP/ds_results.tsv" \
    | awk -v c="$NAME" 'BEGIN{OFS="\t"}{print c,$0}' >> "$COMBINED"

  # mirror the advisory coverage into the combined table (same key columns)
  [[ -s "$CMP/pairs.cover.tsv" ]] \
    && awk -v c="$NAME" 'BEGIN{OFS="\t"}{print c,$0}' "$CMP/pairs.cover.tsv" >> "$COVALL"
done

# fusion breakpoint, as estimate_breakpoint() in figs/fig3/fig3c.R: the juncea_X position minimising the share
# of Anax 12-anchored genes right of it + Anax X-anchored genes left of it; checks, not overwrites, breakpoint.txt
BP_TOL=1000                            # bp; disagreement above this is an error

estimate_breakpoint() {
  local C="$WORKDIR/chrom"
  local GFF="$C/juncea_X.longest.gff3" IDS="$C/juncea_X.ids"
  local NEOP="$WORKDIR/comp/juncea_X__vs__anax_12/mcscanx.pairs.tsv"
  local ANCP="$WORKDIR/comp/juncea_X__vs__anax_X/mcscanx.pairs.tsv"
  local OUT="$WORKDIR/breakpoint.txt"

  local f
  for f in "$GFF" "$IDS" "$NEOP" "$ANCP"; do
    [[ -s "$f" ]] || { echo "[error] breakpoint: missing or empty $f" >&2; return 1; }
  done

  # X-gene midpoints, one line per distinct gene, for each side
  local TMP="$WORKDIR/.bp.$$"
  mkdir -p "$TMP"
  local side
  for side in neo anc; do
    local SRC; [[ "$side" == neo ]] && SRC="$NEOP" || SRC="$ANCP"
    awk -F'\t' -v ids="$IDS" -v gff="$GFF" '
      BEGIN {
        while ((getline k < ids) > 0) XID[k]
        while ((getline l < gff) > 0) {
          n = split(l, g, "\t")
          if (n < 9 || g[3] != "gene") continue
          if (!match(g[9], /ID=[^;]+/)) continue
          MID[substr(g[9], RSTART+3, RLENGTH-3)] = (g[4] + g[5]) / 2
        }
      }
      {
        x = ($1 in XID) ? $1 : (($2 in XID) ? $2 : "")
        if (x == "") next
        sub(/\.t[0-9]+$/, "", x)
        if (x in MID) SEEN[x] = MID[x]
      }
      END { for (g in SEEN) print SEEN[g] }
    ' "$SRC" | sort -n > "$TMP/$side"
    [[ -s "$TMP/$side" ]] || { echo "[error] breakpoint: no $side positions" >&2
                               rm -rf "$TMP"; return 1; }
  done

  # candidate cuts = every anchored position; the best one (ties -> middle of the tied range)
  sort -n -u "$TMP/neo" "$TMP/anc" > "$TMP/cuts"

  local BP
  BP=$(awk -v nf="$TMP/neo" -v af="$TMP/anc" '
    BEGIN {
      while ((getline v < nf) > 0) neo[++nn] = v
      while ((getline v < af) > 0) anc[++na] = v
      best = 1e308
    }
    {
      cx = $1; r = 0; l = 0
      for (i = 1; i <= nn; i++) if (neo[i] >  cx) r++
      for (i = 1; i <= na; i++) if (anc[i] <  cx) l++
      e = r / nn + l / na
      if (e < best - 1e-12) { best = e; lo = cx; hi = cx }
      else if (e <= best + 1e-12)      { hi = cx }
    }
    END { printf "%d\n", int((lo + hi) / 2 + 0.5) }
  ' "$TMP/cuts")

  echo "[breakpoint] estimated $BP bp from $(wc -l < "$TMP/neo") neo and $(wc -l < "$TMP/anc") ancestral anchored X genes"
  rm -rf "$TMP"

  # check against an existing breakpoint.txt (fig3c.R writes it too), else write it
  if [[ -s "$OUT" ]]; then
    local OLD; OLD=$(tr -dc '0-9' < "$OUT")
    local D=$(( BP > OLD ? BP - OLD : OLD - BP ))
    if (( D > BP_TOL )); then
      echo "[error] breakpoint: computed $BP but $OUT holds $OLD (differ by $D bp)." >&2
      echo "[error] One of the two is stale. Delete $OUT to accept the computed value." >&2
      return 1
    fi
    echo "[breakpoint] agrees with existing $OUT ($OLD bp, differ by $D)"
  else
    echo "$BP" > "$OUT"
    echo "[breakpoint] wrote $OUT"
  fi

  return 0
}

estimate_breakpoint \
  || { echo "[error] breakpoint estimation failed" >&2; exit 1; }

# neo-Y -> Anax chr12 miniprot probe of every neo-Y gene (hit = ancestrally present, no hit = candidate
# neo-Y gain); ran in the v8 run, nothing in this repository reads its output
build_entity juncea_Y      # instant: cached above
build_entity anax_12       # instant: cached above
mkdir -p fate/neoY_to_anax12
NOPAIRS="fate/neoY_to_anax12/_probe_all.tsv"; : > "$NOPAIRS"
probe_fate neoY_to_anax12 \
  "$C/juncea_Y.pep.faa" "$C/juncea_Y.longest.gff3" "$C/anax_12.fna" \
  "$NOPAIRS" 1 "${SCAFFOLD[juncea_Y]}" 0 all
echo "[done] neo-Y ancestral-presence probe -> fate/neoY_to_anax12/probe.gff"
