#!/bin/bash
# Rename sequences (make_lastz_chains allows no dots), build 2bit and chrom.sizes, check softmasking
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; activate

# ">OZ369554.1 Aeshna juncea ..." -> ">OZ369554"; <name>.rename.tsv maps old -> new names
for pair in "$REF_NAME:$REF_FA" "$QRY_NAME:$QRY_FA"; do
  n=${pair%%:*}; f=${pair#*:}
  zcat -f "$f" | sed -E '/^>/{s/^>(\S+)\.[0-9]+.*/>\1/}' > "$GEN/$n.fa"
  dup=$(grep '^>' "$GEN/$n.fa" | sort | uniq -d | head -1)
  [ -z "$dup" ] || { echo "duplicate sequence name after rename in $n: $dup"; exit 1; }
  paste <(zcat -f "$f" | grep '^>' | sed 's/^>//; s/ .*//') \
        <(grep '^>' "$GEN/$n.fa" | sed 's/^>//') > "$GEN/$n.rename.tsv"
  awk -v n="$n" '!/^>/{t+=length($0); u+=gsub(/[acgt]/,""); x+=gsub(/[Nn]/,"")} END{printf "%s: %.1f%% lowercase, %.1f%% N\n", n, 100*u/t, 100*x/t; exit (u/t < 0.01)}' "$GEN/$n.fa"
  faToTwoBit "$GEN/$n.fa" "$GEN/$n.2bit"
  twoBitInfo "$GEN/$n.2bit" "$GEN/$n.chrom.sizes"
done

# reference annotation with the same seqid rename (keeps comment lines)
gffread "$REF_GFF" -T -o - | awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} {sub(/\.[0-9]+$/, "", $1); print}' > "$REFIN/$REF_NAME.gtf"
echo "genomes in $GEN; reference GTF in $REFIN"