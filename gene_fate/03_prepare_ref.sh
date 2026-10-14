#!/bin/bash
# Reference annotation -> TOGA2 inputs (BED12, isoforms, U12 introns, CESAR2 profiles)
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; activate
toga2.py prepare-input \
  --ref_2bit "$GEN/$REF_NAME.2bit" \
  --ref_annot "$REFIN/$REF_NAME.gtf" -f gtf \
  --intronic_cores 4 \
  -o "$REFIN"
ls "$REFIN"      # expect: toga.transcripts.bed  toga.isoforms.tsv  toga.U12introns.bed  CESAR2.0_profiles/
