#!/usr/bin/env Rscript
# Collect maps/ into simmaps.rds (read by events_adapter.R) and report transitions per map

suppressPackageStartupMessages({ library(data.table) })
source("config.R")

f <- sort(list.files("maps", pattern = "^simmap_part_[0-9]+_[0-9]+\\.rds$",
                     full.names = TRUE))
if (!length(f)) stop("maps/ is empty -- run the 09c array first", call. = FALSE)

# read every part file (one map or a list of maps), keep only elements that are maps
maps <- list()
for (p in f) {
  x <- readRDS(p)
  maps <- c(maps, if (is.list(x) && is.null(x$maps)) x else list(x))
}
is_map <- logical(length(maps))
for (k in seq_along(maps)) is_map[k] <- !is.null(maps[[k]]$maps)
maps <- maps[which(is_map)]
saveRDS(maps, "simmaps.rds")
cat(sprintf("%d part files -> %d maps in simmaps.rds\n", length(f), length(maps)))

# count transitions per map as events_adapter.R classifies them: count down = fusion
# (SA if the sex class changes), count up = fission, same count + sex-class change = reversal
tab <- matrix(0L, length(maps), 4, dimnames = list(NULL, c("SA", "AA", "fission", "reversal")))
for (i in seq_along(maps)) {
  m <- maps[[i]]
  sa <- aa <- fis <- rv <- 0L
  for (e in seq_along(m$maps)) {
    g <- names(m$maps[[e]]); if (length(g) < 2) next
    for (k in seq_len(length(g) - 1L)) {
      a <- strsplit(g[k], "_", fixed = TRUE)[[1]]
      b <- strsplit(g[k + 1L], "_", fixed = TRUE)[[1]]
      na <- suppressWarnings(as.numeric(a[1])); nb <- suppressWarnings(as.numeric(b[1]))
      if (is.na(na) || is.na(nb)) next
      if (nb < na) { if (a[2] != b[2]) sa <- sa + 1L else aa <- aa + 1L }
      else if (nb > na) fis <- fis + 1L
      else if (a[2] != b[2]) rv <- rv + 1L
    }
  }
  tab[i, ] <- c(sa, aa, fis, rv)
}

# summary across maps
n_ev <- rowSums(tab)
med <- apply(tab, 2, median); rng <- apply(tab, 2, range)
cat(sprintf("transitions per map: median %.0f (range %d-%d)\n",
            median(n_ev), min(n_ev), max(n_ev)))
cat("  broken down as the test counts them:\n")
for (k in colnames(tab))
  cat(sprintf("    %-9s median %5.0f  (range %3d-%3d)\n", k, med[k], rng[1, k], rng[2, k]))
fus <- tab[, "SA"] + tab[, "AA"]
cat(sprintf("  FUSIONS only (SA + AA): median %.0f; SA fraction %.3f\n",
            median(fus), median(tab[, "SA"] / pmax(fus, 1))))
cat(sprintf("  null expectation at n = median: %.3f\n", p_SA(13, "XO")))
cat("  Only SA and AA enter the test. Fissions and reversals are not fusions\n")
cat("  and classify_transition() returns NULL for both.\n")
if (median(tab[, "SA"]) == 0)
  cat("\nWARNING: no SA fusions mapped. Check the sexcls encoding and the root\n",
      "prior.\n", sep = "")
