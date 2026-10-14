#!/usr/bin/env bash
# Align each <marker>.fna (MAFFT), trim to columns present in >= half the taxa (trimAl -gt 0.5),
# and report occupancy
set -eoux pipefail

THREADS="${SLURM_CPUS_PER_TASK:-8}"
mkdir -p align

# stop early if 05_select.R has not written any <marker>.fna
shopt -s nullglob
inputs=()
for f in *.fna; do
  case "$f" in unfiltered.COI.fna|screened.COI.fna) continue ;; esac
  inputs+=("$f")
done
if [ ${#inputs[@]} -eq 0 ]; then
  echo "ERROR: no <marker>.fna files here." >&2
  echo "  06_align_trim.sh aligns the output of 05_select.R. Run, in order:" >&2
  echo "    Rscript 03_screen.R" >&2
  echo "    sbatch  04_blast.sbatch     # wait for hits_<marker>.tsv" >&2
  echo "    Rscript 05_select.R         # writes COI.fna, 16S.fna, ..." >&2
  exit 1
fi
echo "aligning: ${inputs[*]}"

# align (reverse strands adjusted, MAFFT's _R_ prefix stripped), then trim
for f in *.fna; do
  case "$f" in unfiltered.COI.fna|screened.COI.fna) continue ;; esac
  m="${f%.fna}"

  mafft --localpair --maxiterate 1000 --adjustdirectionaccurately \
        --thread "$THREADS" "$f" > "align/${m}.aln.fasta"
  sed -i 's/^>_R_/>/' "align/${m}.aln.fasta"

  trimal -in "align/${m}.aln.fasta" -out "align/${m}.trim.fasta" \
         -gt 0.5 -fasta -htmlout "align/${m}.html"
done

# taxa, columns and occupancy per trimmed alignment -> align/occupancy.tsv
Rscript - <<'RS'
suppressPackageStartupMessages({ library(Biostrings); library(data.table) })
f <- list.files("align", "\\.trim\\.fasta$", full.names = TRUE)
if (!length(f)) {
  cat("no trimmed alignments in align/ -- nothing to summarise\n")
  quit(save = "no")
}
r <- rbindlist(lapply(f, function(p) {
  s <- readBStringSet(p)
  if (!length(s)) return(NULL)
  ch <- as.matrix(letterFrequency(s, "-"))
  data.table(marker = sub("\\.trim\\.fasta$", "", basename(p)),
             taxa = length(s), cols = as.integer(width(s)[1]),
             occupancy = round(1 - sum(ch) / (length(s) * width(s)[1]), 3))
}))
if (!nrow(r)) {
  cat("every trimmed alignment was empty -- check align/*.html from trimAl\n")
  quit(save = "no")
}
print(r[order(-taxa)])
fwrite(r, "align/occupancy.tsv", sep = "\t")
cat("\nOccupancy below ~0.6 means the marker is mostly gaps after trimming and\n")
cat("is adding noise, not signal. Drop it, or raise REF_QUANTILE for it in\n")
cat("config.R and rerun 04/05 with a longer reference.\n")
RS
