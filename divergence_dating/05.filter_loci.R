# Topology / translation / length filters on the 04.align_trees.sbatch output -> filter_report.tsv, locus_table_final.tsv
suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)   # BiocManager::install("Biostrings") if missing
  library(ape)
})

locus_table_path <- "locus_table.tsv"
outdir           <- "per_locus"

min_codons <- 100   # drop loci with post-alignment length below this

# loci with both an alignment and a gene tree
cat("\n=== Stage 3: applying filters ===\n")
locus_table <- fread(locus_table_path)
complete_loci <- locus_table$pg_x[
  file.exists(file.path(outdir, locus_table$pg_x, paste0(locus_table$pg_x, ".best.fas"))) &
  file.exists(file.path(outdir, locus_table$pg_x, paste0(locus_table$pg_x, ".treefile")))
]
cat(sprintf("loci with alignment + tree: %d / %d\n", length(complete_loci), nrow(locus_table)))
fr <- data.table(pg_x = complete_loci)

# Filter 1: neo-X and neo-Y sister in the unrooted ML gene tree (NA if a tip is missing)
sister <- logical(nrow(fr))
for (i in seq_len(nrow(fr))) {
  pg <- fr$pg_x[i]
  tr <- read.tree(file.path(outdir, pg, paste0(pg, ".treefile")))
  if (is.rooted(tr)) tr <- unroot(tr)
  ix <- which(tr$tip.label == "Aeshna_juncea_neoX")
  iy <- which(tr$tip.label == "Aeshna_juncea_neoY")
  if (length(ix) != 1 || length(iy) != 1) {
    sister[i] <- NA
  } else {
    sister[i] <- tr$edge[tr$edge[, 2] == ix, 1] == tr$edge[tr$edge[, 2] == iy, 1]
  }
}
fr[, xy_sister := sister]
cat(sprintf("  Filter 1 (neoX + neoY sister):    %d / %d pass\n",
            sum(fr$xy_sister, na.rm = TRUE), nrow(fr)))

# Filter 2: neo-Y CDS is clean (no frameshift, no internal stop; one terminal stop allowed)
clean  <- logical(nrow(fr))
reason <- character(nrow(fr))
for (i in seq_len(nrow(fr))) {
  pg <- fr$pg_x[i]
  aln <- readDNAStringSet(file.path(outdir, pg, paste0(pg, ".best.fas")))
  y <- aln[names(aln) == "Aeshna_juncea_neoY"]
  if (length(y) == 0) {
    reason[i] <- "neoY not in alignment"
    next
  }
  y_ungapped <- gsub("-", "", as.character(y[[1]]))
  if (nchar(y_ungapped) %% 3 != 0) {
    reason[i] <- "CDS not divisible by 3 (frameshift)"
    next
  }
  prot <- tryCatch(
    as.character(translate(DNAString(y_ungapped), if.fuzzy.codon = "solve")),
    error = function(e) NA_character_
  )
  if (is.na(prot)) {
    reason[i] <- "translation failed"
    next
  }
  if (grepl("\\*", substr(prot, 1, nchar(prot) - 1))) {
    reason[i] <- "internal stop"
    next
  }
  clean[i] <- TRUE
}
fr[, neoY_clean  := clean]
fr[, neoY_reason := reason]
cat(sprintf("  Filter 2 (neoY clean CDS):         %d / %d pass\n",
            sum(fr$neoY_clean), nrow(fr)))

# Filter 3: alignment length >= min_codons
bp <- integer(nrow(fr))
for (i in seq_len(nrow(fr))) {
  pg <- fr$pg_x[i]
  aln <- readDNAStringSet(file.path(outdir, pg, paste0(pg, ".best.fas")))
  bp[i] <- as.integer(width(aln)[1])
}
fr[, aln_bp := bp]
fr[, len_ok := aln_bp >= min_codons * 3]
cat(sprintf("  Filter 3 (aln >= %d codons):     %d / %d pass\n",
            min_codons, sum(fr$len_ok), nrow(fr)))

# keep loci passing all three
fr[, keep := xy_sister & neoY_clean & len_ok]
fr[is.na(keep), keep := FALSE]

# summary
final_loci <- fr[keep == TRUE, pg_x]
kept_bp    <- fr[keep == TRUE, sum(aln_bp)]

cat("\n=== summary ===\n")
cat(sprintf("Final locus set:              %d / %d\n", length(final_loci), nrow(fr)))
cat(sprintf("Total concatenated length:    %s bp\n", format(kept_bp, big.mark = ",")))
cat(sprintf("  = %s codons per position\n", format(kept_bp %/% 3, big.mark = ",")))

if (length(final_loci) < 15)
  cat("WARNING: fewer than 15 loci retained; MCMCtree may be under-powered.\n")
if (kept_bp %/% 3 < 5000)
  cat("WARNING: <5,000 codons per position; deep node ages may be poorly resolved.\n")

# alignment length distribution
cat("\nKept-locus alignment lengths (codons):\n")
kept_codons <- fr[keep == TRUE, aln_bp %/% 3]
print(summary(kept_codons))

# dropped loci and why
if (any(!fr$keep)) {
  cat("\nDropped loci (reasons):\n")
  dropped <- fr[keep == FALSE, .(pg_x, xy_sister, neoY_clean, neoY_reason,
                                 aln_codons = aln_bp %/% 3, len_ok)]
  print(dropped)
}

# outputs
fwrite(fr, "filter_report.tsv", sep = "\t")
locus_table_final <- locus_table[pg_x %in% final_loci]
fwrite(locus_table_final, "locus_table_final.tsv", sep = "\t")
cat("\nWrote: filter_report.tsv, locus_table_final.tsv\n")
cat("Ready for phylip assembly and MCMCtree.\n")
