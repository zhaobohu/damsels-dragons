# Codon-partitioned MCMCtree phylip files from the per-locus PRANK alignments: all loci, then the top dS_XY subsets
suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)
})

locus_table_final_path <- "locus_table_final.tsv"
outdir                 <- "per_locus"
output_full            <- "ajun_concat.phy"
top_fractions          <- c(0.75, 0.5, 0.25)   # dS_XY subsets to emit

# tip order — must match your calibration tree (ajun_calib.tre) exactly
tip_order <- c(
  "Aeshna_juncea_neoX",
  "Aeshna_juncea_neoY",
  "Anax_junius",
  "Somatochlora_arctica",
  "Cordulegaster_boltonii",
  "Sympetrum_striolatum",
  "Orthetrum_coerulescens",
  "Pantala_flavescens",
  "Ischnura_elegans",
  "Platycnemis_pennipes",
  "Archilestes_grandis"
)

# concatenate a set of loci by codon position and write MCMCtree phylip (also used by 09.random_controls.R)
build_phylip <- function(lt, output_phylip, label) {
  cat(sprintf("\n=== Building %s (%d loci) ===\n", label, nrow(lt)))

  concat_by_pos <- list()
  for (pos in c("pos1", "pos2", "pos3")) {
    concat_by_pos[[pos]] <- setNames(rep("", length(tip_order)), tip_order)
  }

  for (i in seq_len(nrow(lt))) {
    pg      <- lt$pg_x[i]
    alnfile <- file.path(outdir, pg, paste0(pg, ".best.fas"))
    aln     <- readDNAStringSet(alnfile)

    # sanity: exactly the expected tips
    if (!setequal(names(aln), tip_order)) {
      missing_tips <- setdiff(tip_order, names(aln))
      extra_tips   <- setdiff(names(aln), tip_order)
      stop(sprintf("locus %s: missing (%s), extra (%s)",
                   pg,
                   paste(missing_tips, collapse = ","),
                   paste(extra_tips,   collapse = ",")))
    }

    # sanity: all taxa same width, codon-multiple
    w <- unique(width(aln))
    if (length(w) != 1) stop(sprintf("locus %s: taxa have differing alignment lengths", pg))
    if (w %% 3 != 0)   stop(sprintf("locus %s: alignment length %d not divisible by 3", pg, w))

    # append codon positions 1, 2, 3 of each tip to the three concatenations
    aln_char <- as.character(aln[tip_order])
    for (tip in tip_order) {
      chars <- strsplit(aln_char[tip], "")[[1]]
      for (p in 1:3) {
        concat_by_pos[[p]][tip] <- paste0(concat_by_pos[[p]][tip],
                                          paste(chars[seq(p, w, by = 3)], collapse = ""))
      }
    }

    if (i %% 5 == 0 || i == nrow(lt)) {
      cat(sprintf("  processed %2d/%2d loci  |  running length = %d bp per position\n",
                  i, nrow(lt), nchar(concat_by_pos$pos1[1])))
    }
  }

  # sanity: equal final lengths across positions and taxa
  for (p in 1:3) {
    lens <- nchar(concat_by_pos[[p]])
    if (length(unique(lens)) != 1)
      stop(sprintf("position %d: taxa have differing final lengths (%s)",
                   p, paste(lens, collapse = ",")))
  }
  L1 <- unname(nchar(concat_by_pos$pos1[1]))
  L2 <- unname(nchar(concat_by_pos$pos2[1]))
  L3 <- unname(nchar(concat_by_pos$pos3[1]))
  if (L1 != L2 || L1 != L3)
    stop(sprintf("codon positions have different lengths: %d / %d / %d", L1, L2, L3))
  L <- L1

  cat(sprintf("  concatenated length per codon position: %s bp\n",
              format(L, big.mark = ",")))
  cat(sprintf("  total across 3 partitions:              %s bp\n",
              format(L * 3, big.mark = ",")))

  # write phylip: 3 sequential blocks, one per codon position
  con <- file(output_phylip, "w")
  for (p in 1:3) {
    cat(sprintf(" %d  %d\n", length(tip_order), L), file = con)
    for (tip in tip_order) {
      cat(sprintf("%s  %s\n", tip, toupper(concat_by_pos[[p]][tip])), file = con)
    }
    cat("\n", file = con)
  }
  close(con)

  cat(sprintf("  wrote: %s\n", output_phylip))
  invisible(L)
}

# locus table (needs dS for the subsets)
lt <- fread(locus_table_final_path)
cat(sprintf("Loaded locus table: %d loci\n", nrow(lt)))
cat(sprintf("Taxa (in phylip order):\n%s\n",
            paste(sprintf("  %s", tip_order), collapse = "\n")))

if (!"dS" %in% names(lt))
  stop("locus_table_final.tsv has no `dS` column - cannot compute top-50% subset")

# BUILD FULL PHYLIP (marker: 09.random_controls.R sources everything up to two lines above it)
build_phylip(lt, output_full, label = "full dataset")

# dS_XY subsets: the highest-dS loci are the most likely to have escaped gene conversion
lt_sorted <- lt[order(-dS)]
outputs   <- output_full

for (tf in top_fractions) {
  pct    <- round(tf * 100)
  n_top  <- max(1, ceiling(nrow(lt_sorted) * tf))
  lt_top <- lt_sorted[seq_len(n_top)]
  outf   <- sprintf("ajun_concat_top%d.phy", pct)

  cat(sprintf("\nTop-%d%% subset: %d loci (dS_XY range %.4f-%.4f)\n",
              pct, n_top, min(lt_top$dS), max(lt_top$dS)))

  build_phylip(lt_top, outf, label = sprintf("top-%d%% dS_XY subset", pct))
  outputs <- c(outputs, outf)
}

# preview the first lines of each alignment
for (path in outputs) {
  cat(sprintf("\n=== preview: %s ===\n", path))
  con <- file(path, "r")
  lines <- readLines(con, n = length(tip_order) + 2)
  close(con)
  for (ln in lines) {
    cat(substr(ln, 1, 80), if (nchar(ln) > 80) " ..." else "", "\n", sep = "")
  }
}
cat(sprintf("\nWrote %d alignments: %s\n",
            length(outputs), paste(outputs, collapse = ", ")))
