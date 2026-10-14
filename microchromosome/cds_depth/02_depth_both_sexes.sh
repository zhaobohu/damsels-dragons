#!/usr/bin/env bash
# Per-CDS read depth for male and female, on the same reference (the male hap1
# assembly), so the two are directly comparable gene by gene.
#
# Usage:
#   ./02_depth_both_sexes.sh \
#       --male-bam  /path/male.hap1.bam \
#       --female-bam /path/female.hap1.bam \
#       --cds-bed   ../../gene_fate/depth_check/cds.bed \
#       [--mapq 20] [--threads 4] [--outdir .] [--also-q0] [--windows 200000]
#
# --also-q0 reruns at MAPQ 0. The gap between the two is itself the repeat
# diagnostic: chr12 gained 13.8% at q0 in the earlier run, the most of any
# chromosome, which is what pointed at mismapping in the first place.

set -euo pipefail

MAPQ=20
THREADS=4
OUTDIR=.
ALSO_Q0=0
WINDOWS=

# parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --male-bam)   MALE="$2"; shift 2 ;;
    --female-bam) FEMALE="$2"; shift 2 ;;
    --cds-bed)    CDS="$2"; shift 2 ;;
    --mapq)       MAPQ="$2"; shift 2 ;;
    --threads)    THREADS="$2"; shift 2 ;;
    --outdir)     OUTDIR="$2"; shift 2 ;;
    --also-q0)    ALSO_Q0=1; shift ;;
    --windows)    WINDOWS="$2"; shift 2 ;;
    -h|--help)    sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

# check required arguments, inputs and mosdepth
: "${MALE:?--male-bam is required}"
: "${FEMALE:?--female-bam is required}"
: "${CDS:?--cds-bed is required}"

for f in "$MALE" "$FEMALE" "$CDS"; do
  [[ -e "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
command -v mosdepth >/dev/null || { echo "mosdepth not on PATH" >&2; exit 1; }

mkdir -p "$OUTDIR"

# cds.bed must be coordinate-sorted for mosdepth --by
SORTED="$OUTDIR/cds.sorted.bed"
sort -k1,1 -k2,2n "$CDS" > "$SORTED"

# per-CDS mean depth with mosdepth
run_one () {  # $1 label  $2 bam  $3 mapq
  local pref="$OUTDIR/$1.cds.q$3"
  echo "[mosdepth] $1 MAPQ$3 -> ${pref}.regions.bed.gz"
  mosdepth -t "$THREADS" -n --mapq "$3" --by "$SORTED" "$pref" "$2"
}

# both sexes at MAPQ, and optionally again at MAPQ 0
run_one male   "$MALE"   "$MAPQ"
run_one female "$FEMALE" "$MAPQ"

if [[ "$ALSO_Q0" -eq 1 && "$MAPQ" -ne 0 ]]; then
  run_one male   "$MALE"   0
  run_one female "$FEMALE" 0
fi

# optional: fixed windows, then whole-chromosome median depth per sex
if [[ -n "$WINDOWS" ]]; then
  for s in male female; do
    bam=$([[ $s == male ]] && echo "$MALE" || echo "$FEMALE")
    echo "[mosdepth] $s windows ${WINDOWS}bp"
    mosdepth -t "$THREADS" -n --mapq "$MAPQ" --by "$WINDOWS" "$OUTDIR/$s.win.q$MAPQ" "$bam"
  done
  echo
  echo "whole-chromosome median depth (windows with depth > 0):"
  for s in male female; do
    zcat "$OUTDIR/$s.win.q$MAPQ.regions.bed.gz" \
    | awk -v s="$s" '$4>0 {v[$1]=v[$1]" "$4; n[$1]++}
        END{for(c in n) print s"\t"c"\t"n[c]"\t"v[c]}' \
    | while read -r sex chr n rest; do
        med=$(printf '%s\n' $rest | sort -g | awk -v n="$n" '{a[NR]=$1} END{print a[int(NR/2)]}')
        printf "%s\t%s\t%s\t%s\n" "$sex" "$chr" "$n" "$med"
      done
  done | sort -k2,2 -k1,1
fi

echo
echo "done. feed the .cds.q${MAPQ}.regions.bed.gz files to 03_mchrom_ratio.R"
