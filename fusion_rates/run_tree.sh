#!/usr/bin/env bash
# Build the tree stage by stage (0-12), checking each stage's output; run inside tmux.
#   bash run_tree.sh [from_stage] [to_stage]    # e.g. "5" resumes from 5, "5 5" runs only 5
set -euo pipefail

# stage range and helpers: run a stage, stop with a resume hint, check outputs, block on sbatch
FROM="${1:-0}"; TO="${2:-12}"
STAGE=0
run() {   # run <n> <label> <command...>
  local n="$1"; shift; local label="$1"; shift
  STAGE="$n"
  if [ "$n" -lt "$FROM" ] || [ "$n" -gt "$TO" ]; then return 0; fi
  echo; echo "================ stage $n: $label ================"; echo
  "$@"
}
die() { echo; echo "STOPPED at stage $STAGE: $*" >&2; echo "resume with:  bash run_tree.sh $STAGE" >&2; exit 1; }
need() { for f in "$@"; do [ -s "$f" ] || die "expected output missing or empty: $f"; done; }

batch() {   # batch <sbatch file> [extra sbatch args...]  -- submit and block
  local sb="$1"; shift
  echo "submitting $sb (blocking until it ends)"
  sbatch --wait "$@" "$sb"
}

trap 'echo; echo "(driver exited at stage $STAGE)"' EXIT

# 0: GBIF name resolution (skipped if name_map.csv is present)
stage0() {
  if [ -s name_map.csv ]; then
    echo "name_map.csv present -- skipping GBIF resolution (delete it to redo)"
  else
    Rscript 00_resolve_names.R
  fi
  need name_map.csv
}
run 0 "resolve names" stage0

# 1: karyotypes, BOLD, backbone -> taxa.tsv
stage1() {
  Rscript 01_prep.R
  need taxa.tsv prep.csv unfiltered.COI.fna backbone_families.tsv
  echo "taxa: $(($(wc -l < taxa.tsv) - 1))   neo-sex: $(awk -F'\t' 'NR>1 && $5!="X"' taxa.tsv | wc -l)"
}
run 1 "prep" stage1

# 2: NCBI gap fill (cached)
stage2() {
  [ -s ncbi_cache.csv ] && echo "ncbi_cache.csv present: $(($(wc -l < ncbi_cache.csv) - 1)) cached sequences"
  Rscript 02_fetch.R
  need prep.csv
  ls raw_*.fasta >/dev/null 2>&1 || die "02_fetch.R wrote no raw_<marker>.fasta"
}
run 2 "fetch" stage2

# 3: COI screen
stage3() {
  Rscript 03_screen.R
  need screened.COI.fna screen_log.tsv
  echo "screened COI records: $(grep -c '^>' screened.COI.fna)"
}
run 3 "screen" stage3

# 4: per-marker reference and BLAST coverage [batch]
stage4() {
  batch 04_blast.sbatch
  need refs/reference_report.tsv
  ls hits_*.tsv >/dev/null 2>&1 || die "BLAST produced no hits_<marker>.tsv"
  echo "references:"; cut -f1-3 refs/reference_report.tsv | column -t | sed 's/^/  /'
  if grep -q "^COI" refs/reference_report.tsv; then
    grep -q "^COI.*folmer" refs/reference_report.tsv \
      || echo "  WARNING: COI reference is NOT the Folmer region -- check MITO_ASM in config.R"
  fi
}
run 4 "blast" stage4

# 5: coverage filter -> <marker>.fna
stage5() {
  Rscript 05_select.R | tee logs/05_select.out
  need tree_taxa.txt coverage_report.tsv
  ls *.fna | grep -v -E '^(unfiltered|screened)\.COI\.fna$' >/dev/null || die "no <marker>.fna written"
}
run 5 "select" stage5

# 6: restrict to MARKERS_TREE
stage6() {
  Rscript 05b_subset.R | tee logs/05b_subset.out
  need tree_taxa.txt
  n=$(wc -l < tree_taxa.txt)
  echo "tips going into the tree: $n"
  [ "$n" -ge 20 ] || die "only $n tips after subsetting -- check MARKERS_TREE"
  grep -q "WARNING: fewer than 5 neo-sex" logs/05b_subset.out && die "focal clade lost its neo-sex tips"
  true
}
run 6 "subset" stage6

# 7: MAFFT + trimAl
stage7() {
  bash 06_align_trim.sh 2>&1 | tail -20
  need align/occupancy.tsv
  ls align/*.trim.fasta >/dev/null 2>&1 || die "no trimmed alignments"
  awk -F'\t' 'NR>1 && $4 < 0.6 {print "  WARNING: " $1 " occupancy " $4 " -- mostly gaps, consider dropping"}' align/occupancy.tsv
}
run 7 "align" stage7

# 8: backbone constraint + partitions
stage8() {
  Rscript 07_constraint.R
  need tree/constraint.tre tree/partitions.nex
}
run 8 "constraint" stage8

# 9: pre-flight check (must exit 0)
stage9() {
  Rscript 07b_check.R | tee logs/07b_check.out || die "07b_check.R found problems -- read logs/07b_check.out"
}
run 9 "pre-flight" stage9

# 10: 10-seed IQ-TREE array [batch]; record the best tree in BEST_TREE
stage10() {
  batch 08_tree.sbatch
  ls run.*/run.*.treefile >/dev/null 2>&1 || die "no treefiles from the IQ-TREE array"
  echo "best scores:"
  grep -H "BEST SCORE" run.*/run.*.log | sort -t: -k3 -g | sed 's/^/  /' | tail -3
  BEST=$(grep -H "BEST SCORE" run.*/run.*.log | sort -t: -k3 -g | tail -1 | cut -d: -f1 | sed 's/\.log$/.treefile/')
  echo "$BEST" > BEST_TREE
  echo "best tree: $BEST"
}
run 10 "tree search" stage10

# 11: UFBoot on the best tree, --wbtl [batch]
stage11() {
  need BEST_TREE
  BEST_TREE="$(cat BEST_TREE)" batch 08b_support.sbatch --export=ALL,BEST_TREE="$(cat BEST_TREE)"
  need support/support.treefile support/support.ufboot
  nboot=$(grep -c ';' support/support.ufboot)
  echo "bootstrap trees: $nboot"
  [ "$nboot" -ge 1000 ] || die "support.ufboot has $nboot trees -- 08b needs --wbtl"
}
run 11 "support" stage11

# 12: chronos on the bootstrap trees [batch], then a root-age check
stage12() {
  rm -f tree/boot_dated_ckpt.rds
  batch 09_date.sbatch
  need tree/odon.dated.tre tree/odon.boot.dated.tre
  Rscript -e '
    suppressPackageStartupMessages(library(ape))
    t <- read.tree("tree/odon.boot.dated.tre"); if (inherits(t,"phylo")) t <- list(t)
    r <- sapply(t, function(x) max(branching.times(x)))
    cat(sprintf("dated: %d trees, root age %.1f (range %.1f-%.1f)\n", length(t), median(r), min(r), max(r)))
    if (median(r) > 400) cat("WARNING: root age > 400 -- the crown calibration did not anchor. Check 09_date.R.\n")'
}
run 12 "date" stage12

# what was built and what comes next
echo
echo "================ tree build complete ================"
echo "  tree/odon.dated.tre         ML tree, dated"
echo "  tree/odon.boot.dated.tre    $(grep -c ';' tree/odon.boot.dated.tre 2>/dev/null || echo '?') dated bootstrap trees"
echo
echo "next: stochastic mapping and the test --"
echo "  Rscript 09b_prep_simmap.R"
echo "  Rscript 09d_make_todo.R <n_total>"
echo "  sbatch --time=0-04:00:00 --array=1-\$(wc -l < todo.txt)%200 09c_simmap.sbatch"
