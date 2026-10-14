#!/usr/bin/env Rscript
# Figures and tables from the stochastic maps: focal z among all candidate clades (A), obs/exp vs
# clade size (B), SA fraction at three levels (C, input to figs/fig2/fig2c.R) and by branch age (D)
#   Rscript 13_figures.R [n_maps]

suppressPackageStartupMessages({ library(ape); library(data.table) })
source("config.R"); source("tree_utils.R"); source("root_tree.R")
source("events_adapter.R")

args  <- commandArgs(trailingOnly = TRUE)
n_use <- if (length(args)) as.integer(args[1]) else min(200L, n_iter_available())

tx <- fread(TAXA_TSV); tx[, tip := gsub(" ", "_", species)]
focal_tips <- gsub(" ", "_", FOCAL_TIPS)

# rows of tr$edge whose child lies below node nd (anc: ancestors_by_node, chld: tr$edge[, 2])
edges_below <- function(nd, anc, chld) {
  below <- logical(length(chld))
  for (k in seq_along(chld)) below[k] <- nd %in% anc[[chld[k]]]
  which(below)
}

rows <- list(); permap <- list(); ages <- list()
for (i in seq_len(n_use)) {
  # tree and fusion events of map i; skip maps with no events or no root
  d <- events_for_iter(i); tr <- d$tree; ev <- as.data.table(d$events)
  if (!is.null(tr) && !is.rooted(tr)) tr <- suppressWarnings(root_tree(tr, tx))
  if (!nrow(ev) || !is.rooted(tr)) next
  ntip <- Ntip(tr); dt <- tips_by_node(tr); anc <- ancestors_by_node(tr)

  # candidate clades: N_MIN..N_MAX tips and no tips shared with the focal clade
  nodes <- which(lengths(dt) >= N_MIN & lengths(dt) <= N_MAX & seq_along(dt) > ntip)
  ft  <- intersect(focal_tips, tr$tip.label)
  foc <- if (length(ft) >= 2) getMRCA(tr, ft) else NA_integer_
  if (!is.na(foc)) {
    fti <- match(ft, tr$tip.label)
    shares <- logical(length(nodes))
    for (k in seq_along(nodes)) shares[k] <- any(dt[[nodes[k]]] %in% fti)
    nodes <- nodes[!shares]
  }

  # skip maps whose focal MRCA exceeds MAX_FOCAL_BLOAT x the focal tip count; add the focal clade
  if (!is.na(foc) && length(dt[[foc]]) > MAX_FOCAL_BLOAT * length(ft)) next
  if (!is.na(foc)) nodes <- union(nodes, foc)
  ev[, p := p_null(nchr, sys)]
  chld <- tr$edge[, 2]

  # per-map SA fractions for panel C (null: every event relabelled SA/AA by its own p_null)
  in_foc <- if (!is.na(foc)) ev$edge %in% edges_below(foc, anc, chld) else rep(FALSE, nrow(ev))
  if (any(ev$edge < 1 | ev$edge > nrow(tr$edge)))
    stop(sprintf("map %d: event edge index outside the tree (%d edges)", i, nrow(tr$edge)), call. = FALSE)
  dep <- node.depth.edgelength(tr); node_age <- max(dep) - dep
  ev[, age := (node_age[tr$edge[edge, 1]] + node_age[tr$edge[edge, 2]]) / 2]
  permap[[length(permap) + 1L]] <- data.table(
    iter = i,
    null_frac = mean(rbinom(nrow(ev), 1L, ev$p)),
    bg_frac   = if (any(!in_foc)) mean(ev$is_SA[!in_foc]) else NA_real_,
    foc_frac  = if (any(in_foc))  mean(ev$is_SA[in_foc])  else NA_real_)
  ages[[length(ages) + 1L]] <- ev[, .(iter = i, age, is_SA, p, in_foc)]

  # z per candidate clade (K_MIN+ fusions, V > 0)
  for (nd in nodes) {
    e <- ev[edge %in% edges_below(nd, anc, chld)]
    if (nrow(e) < K_MIN) next
    E <- sum(e$p); V <- sum(e$p * (1 - e$p))
    if (V == 0) next                    # every fusion in a neo lineage: nothing scorable
    rows[[length(rows) + 1L]] <- data.table(
      iter = i, node = nd, size = length(dt[[nd]]),
      genus = names(sort(table(sub("_.*", "", tr$tip.label[dt[[nd]]])), decreasing = TRUE))[1],
      K = nrow(e), S = sum(e$is_SA), E = E,
      z = (sum(e$is_SA) - E) / sqrt(V), ratio = sum(e$is_SA) / E,
      focal = !is.na(foc) && nd == foc)
  }
}
D <- rbindlist(rows)
if (!nrow(D)) stop("no clades scored -- check K_MIN and the focal clade", call. = FALSE)
fwrite(D, "fig_clade_data.tsv", sep = "\t")
FOC <- D[focal == TRUE]; OTH <- D[focal == FALSE]
lab <- FOCAL_GENUS

# A: histogram of z over the other clades, focal median and 95% range marked
pdf("fig_clade_z.pdf", width = 6.5, height = 4.2)
par(mar = c(4.2, 4.2, 2.2, 1))
h <- hist(OTH$z, breaks = 40, plot = FALSE)
plot(h, col = "grey85", border = "white", main = "",
     xlab = "z  (observed - expected sex-autosome fusions, standardised)",
     ylab = "candidate clades", xlim = range(c(D$z, 0)))
abline(v = 0, col = "grey40", lty = 3)
abline(v = median(FOC$z), col = "#D85A30", lwd = 2.5)
q <- quantile(FOC$z, c(.025, .975))
rect(q[1], 0, q[2], max(h$counts) * 0.06, col = "#D85A3033", border = NA)
text(median(FOC$z), max(h$counts) * 0.95, sprintf("  %s", lab),
     col = "#D85A30", adj = 0, font = 2)
text(median(FOC$z), max(h$counts) * 0.86,
     sprintf("  z = %.2f, rank %d of %d", median(FOC$z),
             round(median(D[, .(r = rank(-z)[focal][1]), by = iter]$r, na.rm = TRUE)),
             round(median(D[, .N, by = iter]$N))), col = "#D85A30", adj = 0, cex = 0.8)
dev.off()

# B: obs/exp vs clade size; other clades one point per (genus, size bin), focal clade once at its
# median with bars for its size and ratio ranges (as.numeric: data.table needs one type per column)
pdf("fig_obsexp_size.pdf", width = 6.5, height = 4.2)
par(mar = c(4.2, 4.2, 2.2, 1))
agg <- OTH[, .(size = as.numeric(median(size)), ratio = as.numeric(median(ratio)),
               z = as.numeric(median(z)), n = .N),
           by = .(genus, node_grp = round(size / 3))]
fa <- FOC[, .(size = as.numeric(median(size)), lo = min(size), hi = max(size),
              ratio = as.numeric(median(ratio)),
              rlo = quantile(ratio, .025), rhi = quantile(ratio, .975))]
# clades with no SA fusions (ratio 0) are drawn as x at a floor below the smallest positive ratio
pos  <- agg$ratio[agg$ratio > 0]
flr  <- if (length(pos)) min(c(pos, fa$rlo[fa$rlo > 0])) / 2 else 0.1
agg[, plot_y := pmax(ratio, flr)]
yl   <- range(c(agg$plot_y, fa$rlo, fa$rhi, 1)); yl[1] <- max(yl[1], flr * 0.8)
plot(agg$size, agg$plot_y, pch = ifelse(agg$ratio == 0, 4, 21),
     bg = "grey80", col = "grey60", cex = 0.9, log = "y",
     xlim = range(c(agg$size, fa$lo, fa$hi)), ylim = yl,
     xlab = "clade size (tips)", ylab = "observed / expected SA fusions", main = "")
abline(h = 1, col = "grey40", lty = 3)
if (any(agg$ratio == 0))
  mtext(sprintf("x = %d clades with no SA fusions", sum(agg$ratio == 0)),
        side = 1, line = -1.2, adj = 0.02, cex = 0.7, col = "grey45")
top <- agg[order(-ratio)][1:min(4, .N)]
text(top$size, top$plot_y, top$genus, cex = 0.7, pos = 4, col = "grey30")
segments(fa$lo, fa$ratio, fa$hi, fa$ratio, col = "#D85A30", lwd = 1.5)
segments(fa$size, fa$rlo, fa$size, fa$rhi, col = "#D85A30", lwd = 1.5)
points(fa$size, fa$ratio, pch = 21, bg = "#D85A30", col = "#D85A30", cex = 1.7)
text(fa$size, fa$ratio, lab, cex = 0.85, pos = 4, col = "#D85A30", font = 2, offset = 0.8)
dev.off()

# C: SA fraction at three levels (null, rest of tree, focal) as distributions across maps
P <- rbindlist(permap)
fwrite(P, "fig_levels_permap.tsv", sep = "\t")
xs <- list(P$null_frac, P$bg_frac[!is.na(P$bg_frac)], P$foc_frac[!is.na(P$foc_frac)])
med <- lo <- hi <- numeric(3)
for (k in 1:3) {
  med[k] <- median(xs[[k]])
  lo[k]  <- quantile(xs[[k]], .025)
  hi[k]  <- quantile(xs[[k]], .975)
}
lv <- data.table(level = c("random\nexpectation", "rest of\nthe tree", lab),
                 median = med, lo = lo, hi = hi)
fwrite(lv, "fig_levels.tsv", sep = "\t")
pdf("fig_levels.pdf", width = 5, height = 4.2)
par(mar = c(4.2, 4.6, 1.5, 1))
cols <- c("#BFBFBF", "#8C8C8C", "#D85A30")   # hex so an alpha suffix is valid
plot(NA, xlim = c(0.5, 3.5), ylim = c(0, max(unlist(xs)) * 1.08), xaxt = "n",
     xlab = "", ylab = "sex-autosome fusions / all fusions")
axis(1, at = 1:3, labels = lv$level, padj = 0.6, tick = FALSE)
for (k in 1:3) {
  x <- xs[[k]]; d <- density(x, adjust = 1.2, from = 0, to = 1)
  w <- 0.42 * d$y / max(d$y)
  polygon(c(k - w, rev(k + w)), c(d$x, rev(d$x)), col = paste0(cols[k], "AA"), border = cols[k])
  segments(k - 0.18, median(x), k + 0.18, median(x), lwd = 2.5, col = "grey15")
}
text(1:3, lv$hi + 0.02, sprintf("%.3f", lv$median), cex = 0.85)
dev.off()

# D: SA fraction vs null by branch age, in equal-count age bins (unique() guards against ties)
A <- rbindlist(ages)
brk <- unique(quantile(A$age, seq(0, 1, 0.125)))
if (length(brk) < 3) brk <- unique(quantile(A$age, c(0, .5, 1)))
A[, bin := cut(age, breaks = brk, include.lowest = TRUE)]
byage <- A[, .(mid = median(age), n = .N, obs = mean(is_SA), null = mean(p)), by = bin][order(mid)]
fwrite(byage, "fig_by_age.tsv", sep = "\t")
pdf("fig_by_age.pdf", width = 6, height = 4.2)
par(mar = c(4.2, 4.6, 1.5, 1))
plot(byage$mid, byage$obs, type = "b", pch = 19, col = "grey20", lwd = 2,
     ylim = c(0, max(byage$obs) * 1.15), xlim = rev(range(byage$mid)),
     xlab = "branch age (Ma, root at left)", ylab = "sex-autosome fusions / all fusions")
lines(byage$mid, byage$null, type = "b", pch = 1, col = "grey60", lwd = 1.5, lty = 2)
legend("topleft", c("mapped fusions", "random expectation"), col = c("grey20", "grey60"),
       lwd = c(2, 1.5), lty = c(1, 2), pch = c(19, 1), bty = "n", cex = 0.85)
dev.off()

# summary (byage is youngest bin first)
cat(sprintf("clades scored: %d across %d maps\n", nrow(D), uniqueN(D$iter)))
cat(sprintf("focal z median %.2f; other clades median %.2f\n", median(FOC$z), median(OTH$z)))
cat(sprintf("SA fraction -- null %.3f, background %.3f, focal %.3f  (medians across maps)\n",
            lv$median[1], lv$median[2], lv$median[3]))
cat("SA fraction by branch-age bin (youngest first, midpoint in Ma):\n")
cat(sprintf("  %8s %8s %8s\n", "age", "mapped", "null"))
for (k in seq_len(nrow(byage)))
  cat(sprintf("  %8.1f %8.3f %8.3f\n", byage$mid[k], byage$obs[k], byage$null[k]))
cat("written: fig_clade_z.pdf, fig_obsexp_size.pdf, fig_levels.pdf, fig_by_age.pdf\n")
cat("         + fig_clade_data.tsv, fig_levels.tsv, fig_levels_permap.tsv, fig_by_age.tsv\n")
