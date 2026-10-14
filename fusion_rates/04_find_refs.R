#!/usr/bin/env Rscript
# One BLAST reference per marker (common window for the coverage filter) -> refs/<marker>.ref.fna
# COI: Folmer region cut from the A. juncea mitogenome; other markers: a short record from the data

suppressPackageStartupMessages({ library(Biostrings); library(data.table) })
source("config.R")

dir.create("refs", showWarnings = FALSE)

# Folmer primers
lco <- DNAString("GGTCAACAAATCATAAAGATATTGG")    # LCO1490
hco <- DNAString("TAAACTTCAGGGTGACCAAAAAATCA")   # HCO2198

# region between the primers if exactly one ~658 bp product is found, else NULL
cut_folmer <- function(s, mm) {
  f <- matchPattern(lco, s, max.mismatch = mm)
  r <- matchPattern(reverseComplement(hco), s, max.mismatch = mm)
  if (!length(f) || !length(r)) return(NULL)
  p <- expand.grid(i = seq_along(f), j = seq_along(r))
  p$st <- end(f)[p$i] + 1L
  p$en <- start(r)[p$j] - 1L
  p <- p[abs(p$en - p$st + 1L - 658L) <= 12L, ]
  if (nrow(p) != 1L) return(NULL)          # ambiguous -- refuse rather than guess
  subseq(s, p$st, p$en)
}

# COI reference: Folmer cut on one strand of the mitogenome, kept if it has a stop-free frame
coi_ref <- NULL
if (file.exists(MITO_ASM)) {
  g <- readDNAStringSet(MITO_ASM)
  names(g) <- sub(" .*", "", names(g))
  if (MITO_SEQ %in% names(g)) {
    s  <- g[[MITO_SEQ]]
    fw <- cut_folmer(s, PRIMER_MM)
    rv <- cut_folmer(reverseComplement(s), PRIMER_MM)
    if (xor(is.null(fw), is.null(rv))) {
      coi_ref <- if (is.null(fw)) rv else fw
      ok <- logical(3)
      for (f in 0:2) {
        w  <- (length(coi_ref) - f) %/% 3L * 3L
        aa <- as.character(translate(subseq(coi_ref, f + 1L, width = w),
                genetic.code = getGeneticCode("5"), if.fuzzy.codon = "solve"))
        ok[f + 1] <- !grepl("\\*", sub("\\*+$", "", aa))
      }
      if (!any(ok)) { warning("Folmer cut has no stop-free frame"); coi_ref <- NULL }
    }
  } else warning("contig ", MITO_SEQ, " not in ", MITO_ASM)
} else warning(MITO_ASM, " not found")

# write one reference per marker
rep <- list()
for (m in MARKERS) {
  out <- sprintf("refs/%s.ref.fna", m)

  # COI from the Folmer cut
  if (m == "COI" && !is.null(coi_ref)) {
    r <- DNAStringSet(coi_ref)
    names(r) <- paste0("COI5P_folmer_", MITO_SEQ)
    writeXStringSet(r, out)
    rep[[length(rep) + 1]] <- data.table(marker = m, method = "folmer",
                                         ref_len = width(r), n_records = NA_integer_)
    next
  }

  # otherwise the non-empty record at length quantile q (REF_QUANTILE or the default)
  f <- if (m == "COI") "screened.COI.fna" else sprintf("raw_%s.fasta", m)
  if (!file.exists(f) || file.size(f) == 0) { cat("skip", m, "(no records)\n"); next }
  q <- if (!is.null(REF_QUANTILE[[m]])) REF_QUANTILE[[m]] else REF_QUANTILE_DEFAULT
  recs <- readDNAStringSet(f)
  recs <- recs[width(recs) > 0]
  if (!length(recs)) next
  r <- recs[order(width(recs))[max(1L, ceiling(q * length(recs)))]]
  names(r) <- sprintf("%s_ref_q%02d_%s", m, as.integer(100 * q),
                      gsub("[^A-Za-z0-9]+", "_", names(r)))
  writeXStringSet(r, out)

  # best possible retention: records shorter than SC_MIN * ref_len cannot pass step 05
  w <- width(readDNAStringSet(f))
  ceiling_pct <- mean(w >= SC_MIN * width(r))
  rep[[length(rep) + 1]] <- data.table(marker = m,
    method = sprintf("data q%.2f", q), ref_len = width(r), n_records = length(w),
    max_retain = round(ceiling_pct, 3))
}

# report, warning about markers that cannot retain 90% of records
rp <- rbindlist(rep, fill = TRUE)
print(rp)
low <- rp[!is.na(max_retain) & max_retain < 0.9]
if (nrow(low)) {
  cat("\nWARNING: these markers cannot retain >=90% of records at SC_MIN =",
      SC_MIN, "\n", sep = " ")
  print(low[, .(marker, ref_len, max_retain)])
  cat("The reference is too long for the shorter records. LOWER REF_QUANTILE\n")
  cat("for these markers in config.R -- do not lower SC_MIN, which would let\n")
  cat("genuinely non-overlapping records through.\n")
}
fwrite(rp, "refs/reference_report.tsv", sep = "\t")
if (is.null(coi_ref))
  cat("\nWARNING: COI fell back to a data-driven reference. The Folmer cut is\n",
      "the principled one -- check MITO_ASM / MITO_SEQ / PRIMER_MM in config.R.\n", sep = "")
