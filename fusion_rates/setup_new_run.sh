#!/usr/bin/env bash
# Create a fresh run directory for a tree build, carrying over the download cache and name fixes.
#   bash setup_new_run.sh <old_run_dir> <new_run_dir>
set -euo pipefail

OLD="${1:?usage: setup_new_run.sh <old_run_dir> <new_run_dir>}"
NEW="${2:?usage: setup_new_run.sh <old_run_dir> <new_run_dir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[ -d "$OLD" ] || { echo "old run dir not found: $OLD" >&2; exit 1; }
[ -e "$NEW" ] && { echo "refusing to overwrite existing: $NEW" >&2; exit 1; }
mkdir -p "$NEW/align" "$NEW/tree" "$NEW/logs"

# scripts; the old dir's reconcile.R (your name fixes) wins over the shipped one
cp "$HERE"/*.R "$HERE"/*.sh "$HERE"/*.sbatch "$NEW"/
[ -f "$OLD/reconcile.R" ] && cp "$OLD/reconcile.R" "$NEW/reconcile.R"

# cached downloads. Not name_map.csv (stage 0 must rerun so reconcile.R overrides apply)
# and not raw_*.fasta (02_fetch.R rebuilds them from the cache for the markers in use).
for f in ncbi_cache.csv ncbi_alias.csv; do
  if [ -f "$OLD/$f" ]; then cp "$OLD/$f" "$NEW/"; echo "  copied $f"; fi
done
echo "  NOT copying name_map.csv: stage 0 must rerun so reconcile.R changes apply"

cd "$NEW"
echo "  markers that will be fetched/BLASTed: $(Rscript -e "source('config.R'); cat(paste(MARKERS, collapse=', '))")"

# check that the data paths in config.R resolve
missing=0
for v in KARYO_CSV BOLD_TSV BACKBONE_TRE; do
  p="$(Rscript -e "source('config.R'); cat($v)" 2>/dev/null)"
  if [ -f "$p" ]; then echo "  ok   $v -> $p"; else echo "  MISSING $v -> $p"; missing=1; fi
done
[ "$missing" -eq 0 ] || echo "  fix the paths in config.R before running" >&2

echo
echo "next:  cd $NEW && bash run_tree.sh"
