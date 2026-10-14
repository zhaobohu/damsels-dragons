#!/usr/bin/env Rscript
# Leave-one-out on the focal tips across every tree: does dropping one taxon recover the focal clade?
#   Rscript check_focal_monophyly.R
suppressPackageStartupMessages({ library(ape); library(data.table) })
source("config.R"); source("tree_utils.R"); source("root_tree.R")

tx <- fread(TAXA_TSV)
foc <- gsub(" ", "_", FOCAL_TIPS)

# the run.N search trees and the support tree
rd <- list.dirs(".", recursive = FALSE, full.names = FALSE)
rd <- rd[grepl("^run\\.[0-9]+$", rd)]
tf <- c(sprintf("%s/%s.treefile", rd, rd), "support/support.treefile")
tf <- tf[file.exists(tf)]
if (!length(tf)) stop("no trees found", call. = FALSE)

# per tree: focal MRCA size / focal tips (bloat) with all focal tips and with each one dropped
res <- list()
for (p in tf) {
  tr <- read.tree(p); tr$tip.label <- sub("^_R_", "", tr$tip.label)
  if (!is.rooted(tr)) tr <- suppressWarnings(root_tree(tr, tx))
  ft <- intersect(foc, tr$tip.label)
  if (length(ft) < 3 || !is.rooted(tr)) next
  d <- tips_by_node(tr)
  full <- length(d[[getMRCA(tr, ft)]])
  one <- setNames(numeric(length(ft)), ft)
  for (dropped in ft) {
    r <- setdiff(ft, dropped)
    one[dropped] <- length(d[[getMRCA(tr, r)]]) / length(r)
  }
  res[[length(res) + 1]] <- data.table(tree = basename(p), n_focal = length(ft),
                                       bloat_all = full / length(ft), drop = names(one), bloat_without = one)
}
res <- rbindlist(res)

cat(sprintf("%d trees, %d focal tips\n", uniqueN(res$tree), res$n_focal[1]))
cat(sprintf("bloat with all focal tips: %.2f-%.2f\n\n",
            min(res$bloat_all), max(res$bloat_all)))
s <- res[, .(trees = .N, min = round(min(bloat_without), 2),
             max = round(max(bloat_without), 2),
             n_mono = sum(bloat_without <= 1.001)), by = drop][order(max)]
cat("leave-one-out, across every tree:\n")
print(s, nrows = 40)

# judge on the worst bloat after removal and the gap to the runner-up
best <- s[1]; second <- if (nrow(s) > 1) s[2] else NULL
gap  <- if (!is.null(second)) second$min - best$max else Inf
usable <- best$max <= MAX_FOCAL_BLOAT

cat(sprintf("\nbest single removal: %s -- bloat %.2f-%.2f across %d trees",
            best$drop, best$min, best$max, best$trees))
cat(sprintf(" (clean in %d)\n", best$n_mono))
if (!is.null(second))
  cat(sprintf("next best        : %s -- bloat %.2f-%.2f\n",
              second$drop, second$min, second$max))

if (usable && gap > 0.3) {
  cat("\nDECISIVE: that one taxon accounts for the non-monophyly, and removing it\n")
  cat(sprintf("keeps bloat under MAX_FOCAL_BLOAT (%.2f) in EVERY tree, so no map\n", MAX_FOCAL_BLOAT))
  cat("would be refused. The gap to the next candidate is large, so this is not\n")
  cat("a tie being broken arbitrarily.\n")
  cat("  Still act on it ONLY with independent taxonomic support. Dropping a taxon\n")
  cat("  because it breaks monophyly is circular; dropping one that other molecular\n")
  cat("  studies already place outside the genus is not.\n")
} else if (usable) {
  cat("\nUSABLE but not decisive: removal keeps bloat under MAX_FOCAL_BLOAT, but\n")
  cat("another candidate is nearly as good, so the choice is not forced by the data.\n")
} else {
  cat("\nNo single removal brings bloat under MAX_FOCAL_BLOAT in every tree.\n")
  cat("Removing taxa one at a time until the clade appears would be circular.\n")
}
