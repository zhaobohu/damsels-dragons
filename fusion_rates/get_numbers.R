#!/usr/bin/env Rscript
# Remaining manuscript numbers from files already in the run directory (nothing is rerun or written)
#   Rscript get_numbers.R          # sections 1-3
#   Rscript get_numbers.R clades   # also names the clades in section 4 (loads simmaps.rds)

suppressPackageStartupMessages({ library(ape); library(data.table) })
do_clades <- "clades" %in% commandArgs(trailingOnly = TRUE)
q3 <- function(x, d = 2) sprintf(paste0("%.", d, "f (95%% %.", d, "f-%.", d, "f)"),
                                 median(x), quantile(x, .025), quantile(x, .975))

# 1. fitted rates per rate class, including the reversal rate 09b does not print
pr <- readRDS("prep.rds")
cls <- c("AA fusion", "SA fusion", "fission", "reversal")
cat("1. fitted rates (per lineage per Myr):\n")
for (k in seq_along(cls)) {
  v <- unique(signif(pr$Qfit[pr$tm == k], 8))
  cat(sprintf("     %-10s %s%s\n", cls[k],
              if (length(v)) sprintf("%.4g", v[1]) else "not in model",
              if (length(v) > 1) "  (WARNING: more than one value)" else ""))
}

# 2. Aeshna crown clade across bootstrap trees, and the result split by exactly-15 crowns
P <- fread("timeslice_permap.tsv")[T_factor == 1]
cat(sprintf("2. Aeshna crown clade, %d maps: exactly 15 species in %.1f%%, >30 in %d; crown age %s Ma\n",
            nrow(P), 100 * mean(P$focal_tips == 15), sum(P$focal_tips > 30), q3(P$T_Ma, 1)))
print(table(size = cut(P$focal_tips, c(14, 15, 20, 30, Inf),
                       labels = c("15", "16-20", "21-30", ">30"))))
cat("   split by whether the crown clade is exactly the 15 Aeshna:\n")
print(P[, .(maps = .N, crown_age = round(median(T_Ma), 1), top_z = round(mean(rank_z == 1), 3),
            z = round(median(focal_z), 2), p_maxT = round(median(p_maxT), 4)),
        by = .(only_15 = focal_tips == 15)], row.names = FALSE)

# 3. focal E and obs/exp from the time-slice run
C1 <- fread("timeslice_clades.tsv")[T_factor == 1]
Fc <- C1[focal == TRUE]
cat(sprintf("3. focal clade: K %.0f, S %.0f, E %.2f (medians); S/E %s\n",
            median(Fc$K), median(Fc$S), median(Fc$E), q3(Fc$S / Fc$E)))

# 4. clades with higher z than Aeshna in the same map
fz <- Fc[, .(iter, fz = z)]
B  <- merge(C1[focal == FALSE], fz, by = "iter")[z > fz]
cat(sprintf("4. %d of %d maps have a clade with higher z than Aeshna\n", uniqueN(B$iter), nrow(P)))
if (!do_clades) {
  cat("   rerun with 'clades' to name them\n")
} else if (nrow(B)) {
  maps <- readRDS("simmaps.rds")
  tx <- fread("taxa.tsv"); tx[, tip := gsub(" ", "_", species)]

  # label each clade by its top 3 genera and its families; edge indexes maps[[iter]]$edge
  clade_lab <- character(nrow(B))
  for (r in seq_len(nrow(B))) {
    m <- maps[[B$iter[r]]]; m$tip.label <- sub("^_R_", "", m$tip.label)
    nd <- m$edge[B$edge[r], 2]
    if (nd <= Ntip(m)) {
      sp <- m$tip.label[nd]
    } else {
      # walk down from nd collecting tips
      d <- nd; tip_idx <- integer(0)
      while (length(d)) { ch <- m$edge[m$edge[, 1] %in% d, 2]; tip_idx <- c(tip_idx, ch[ch <= Ntip(m)]); d <- ch[ch > Ntip(m)] }
      sp <- m$tip.label[tip_idx]
    }
    g  <- sort(table(sub("_.*", "", sp)), decreasing = TRUE)
    f  <- unique(na.omit(tx$family[match(sp, tx$tip)]))
    clade_lab[r] <- sprintf("%s%s [%s]", paste(head(names(g), 3), collapse = "/"),
                            if (length(g) > 3) "/..." else "", paste(f, collapse = "/"))
  }
  B[, clade := clade_lab]

  # one row per clade; pct is computed inside the aggregation (data.table skips the first print()
  # after a := on the same table) and medians are as.numeric (integer vs double across groups)
  out <- B[, .(maps = uniqueN(iter), pct_of_maps = round(100 * uniqueN(iter) / nrow(P), 1),
               species = as.numeric(median(n_tips)), fusions = as.numeric(median(K)),
               z = round(median(z), 2)),
           by = clade][order(-maps)]
  print(out, row.names = FALSE)
}
