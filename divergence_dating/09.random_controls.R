#!/usr/bin/env Rscript
# Size-matched random locus subsets (same locus count, total length within tol) as controls for the dS_XY subsets
library(data.table)
library(Biostrings)

setwd("/path/to/project/divergence_dating")
n_rep     <- 999   # permutation p = (1 + k) / (n_rep + 1) comes in steps of 0.001
fractions <- c(0.75, 0.5, 0.25)
tol       <- 0.15
seed      <- 42

# IDs are rand<pct>_<3 digits>; 10.controls.sbatch and 11.collect.sh match them as rand??_???
if (n_rep > 999) stop("n_rep > 999 breaks the rand??_??? IDs used downstream")

# refuse to mix control sets: 11.collect.sh would pick up a leftover rand??_??? run as part of this one
old <- c(Sys.glob("rand??_???"), Sys.glob("ajun_concat_rand*.phy"),
         if (file.exists("controls_manifest.tsv")) "controls_manifest.tsv")
if (length(old))
  stop(length(old), " leftover item(s) from a previous control set, e.g. ", old[1],
       "\n  Move them aside first.")

# build_phylip(), tip_order and lt: 06.build_phylip.R up to its '# BUILD FULL PHYLIP' marker
src <- readLines("06.build_phylip.R")
cut <- grep("^# BUILD FULL PHYLIP", src)
if (!length(cut)) stop("marker '# BUILD FULL PHYLIP' not found in 06.build_phylip.R")
eval(parse(text = paste(src[1:(cut[1] - 2)], collapse = "\n")))
if (!exists("build_phylip") || !exists("lt"))
  stop("sourcing 06.build_phylip.R did not define build_phylip() and lt")

# alignment length of each locus, then loci by decreasing dS
aln_len <- numeric(nrow(lt))
for (i in seq_len(nrow(lt))) {
  pg <- lt$pg_x[i]
  aln_len[i] <- sum(width(readDNAStringSet(file.path("per_locus", pg,
                                                     paste0(pg, ".best.fas"))))[1])
}
lt[, len := aln_len]

ds_col <- if ("dS" %in% names(lt)) "dS" else grep("^dS", names(lt), value = TRUE)[1]
setorderv(lt, ds_col, order = -1L)

man <- data.table(id = character(), frac = numeric(), phy = character(),
                  nloci = integer(), nt = integer(),
                  ds_min = numeric(), ds_max = numeric(), med_len = numeric())

# n_rep draws per fraction with the same locus count and total length within tol
set.seed(seed)
for (tf in fractions) {
  pct    <- round(tf * 100)
  n_top  <- max(1L, ceiling(nrow(lt) * tf))
  target <- sum(lt[seq_len(n_top)]$len)
  cat(sprintf("\n=== top-%d%%: %d loci, %d nt, median locus %.0f nt\n",
              pct, n_top, target, median(lt[seq_len(n_top)]$len)))

  keys <- character(n_rep)             # locus set of each draw, to count repeats
  for (r in seq_len(n_rep)) {
    ok <- FALSE
    for (try in 1:50000) {
      s <- lt[sample(.N, n_top)]
      if (abs(sum(s$len) - target) <= tol * target) { ok <- TRUE; break }
    }
    if (!ok) stop("no size-matched draw for top-", pct, "%, rep ", r)
    keys[r] <- paste(sort(s$pg_x), collapse = ",")
    id   <- sprintf("rand%d_%03d", pct, r)
    outf <- sprintf("ajun_concat_%s.phy", id)
    build_phylip(s, outf, label = id)
    man <- rbind(man, data.table(id = id, frac = tf, phy = outf,
                                 nloci = nrow(s), nt = sum(s$len),
                                 ds_min = min(s[[ds_col]]),
                                 ds_max = max(s[[ds_col]]),
                                 med_len = median(s$len)))
    if (r %% 100 == 0) cat(sprintf("  %d / %d\n", r, n_rep))
  }
  cat(sprintf("  %d controls written, %d distinct locus sets\n",
              n_rep, uniqueN(keys)))
}

fwrite(man, "controls_manifest.tsv", sep = "\t")
cat(sprintf("\nmanifest: %d rows -> controls_manifest.tsv\n", nrow(man)))
