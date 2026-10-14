#!/bin/bash
#SBATCH -J toga2
#SBATCH -p single
#SBATCH -N 1 -n 1 -c 64 --mem-per-cpu=3500M
#SBATCH -t 7-00:00:00
#SBATCH -o toga2.%j.out
# Whole TOGA2 run inside this one job (Nextflow local executor); resume with the same -o plus --resume_from <step>
# steps: setup feature_extraction classification preprocessing aggregate_preprocessing_res alignment
#        aggregate_cesar_res gene_inference loss_summary orthology summarize_trees finalize ucsc_report
set -euo pipefail
source /path/to/project/gene_fate/config.sh; activate
CHAIN=$(ls "$CHAINS"/results/07_final/*.chain.gz | head -1)
echo "chain: $CHAIN"

# one Nextflow config per parallel step; the five CESAR memory bins run at once, so 1 slot each for the
# 5-15 GB bins and the rest of the budget to the 3 GB bin
NFC=$RUN/toga2_nfconfig; mkdir -p "$NFC"
nfcfg() {   # <file> <slots> <GB per task>
  cat > "$NFC/$1" <<CFG
process.executor   = 'local'
process.cpus       = 1
process.memory     = '$3 GB'
process.time       = '48h'
executor.cpus      = $2
executor.memory    = '$(( $2 * $3 )) GB'
executor.queueSize = $2
CFG
}
min() { echo $(( $1 < $2 ? $1 : $2 )); }
nfcfg preprocessing.config $(min $NCPU $(( NMEM / 6 ))) 6      # TOGA2's default 6 GB per task
nfcfg orthology.config     $(min $NCPU $(( NMEM / 6 ))) 6
for b in 5 7 10 15; do nfcfg alignment_$b.config 1 $b; done
b3=$(min $(( NCPU - 4 )) $(( (NMEM - 37) / 3 )))
nfcfg alignment_3.config $b3 3
echo "alignment slots: 3 GB bin = $b3, 5/7/10/15 GB bins = 1 each"

# TOGA2's human non-canonical U2 profiles: prepare-input wrote none (no such Anax introns), but run needs them
HUMAN_PROF=$TOGA2/supply/CESAR2.0/profiles/human
ls "$HUMAN_PROF"/nonCanon_U2_acceptor.tsv "$HUMAN_PROF"/nonCanon_U2_donor.tsv >/dev/null

# --annotate_paralogs aligns paralog-labelled chains too (e.g. a neo-Y copy whose neo-X copy is the ortholog)
toga2.py run \
  --ref_2bit   "$GEN/$REF_NAME.2bit" \
  --query_2bit "$GEN/$QRY_NAME.2bit" \
  --chain_file "$CHAIN" \
  --ref_annotation   "$REFIN/toga.transcripts.bed" \
  --isoform_file     "$REFIN/toga.isoforms.tsv" \
  --u12_file         "$REFIN/toga.U12introns.bed" \
  --cesar_profile_dir "$REFIN/CESAR2.0_profiles" \
  --cesar_non_canon_u2_acceptor "$HUMAN_PROF/nonCanon_U2_acceptor.tsv" \
  --cesar_non_canon_u2_donor    "$HUMAN_PROF/nonCanon_U2_donor.tsv" \
  --no_spliceai \
  --annotate_paralogs \
  --parallel_strategy local --nextflow_config_dir "$NFC" \
  --keep_nextflow_log \
  --project_name "${REF_NAME}_to_${QRY_NAME}" \
  -o "$OUT" -v
cat "$OUT/summary.txt"
