#!/usr/bin/env Rscript
# Fig. 2D: per map, the focal clade's z minus the best competing equal-age clade's z

suppressPackageStartupMessages(library(data.table))

CLADES <- "../../fusion_rates/timeslice_clades.tsv"
PERMAP <- "../../fusion_rates/timeslice_permap.tsv"   # for p_maxT
OUT    <- "d.aeshna.rank.svg"
FOCAL_LABEL <- "Aeshna"
T_USE  <- 1.0

WIDTH <- 30; HEIGHT <- 9
COL_FOC <- "#B22222"; COL_OTH <- "#386cb0"
N_BREAKS <- 45
AXIS_LWD <- 5; AXIS_CEX <- 2.2; AXIS_TCK <- -0.04; AXIS_NUM_LINE <- 2.6
LAB_CEX <- 2.4; LAB_LINE <- 5.2; YLAB_LINE <- 7.5

# per map: focal z, the largest z among the other slice clades, and the margin between them
C <- fread(CLADES)[T_factor == T_USE]
G <- C[, .(fz = z[focal][1],
           best = if (any(!focal)) max(z[!focal], na.rm = TRUE) else NA_real_), by = iter]
G <- G[is.finite(fz) & is.finite(best)]
G[, margin := fz - best]

# histogram of the margin as a density: red where the focal clade is ahead, blue where not
ax <- function(side, ...) axis(side, lwd = AXIS_LWD, cex.axis = AXIS_CEX, tck = AXIS_TCK,
                               col = "grey25", col.axis = "grey15",
                               mgp = c(3, AXIS_NUM_LINE, 0), ...)
svg(OUT, width = WIDTH, height = HEIGHT)
par(mar = c(LAB_LINE + 3, YLAB_LINE + 2.5, 1.5, 1.5), xpd = NA)
h <- hist(G$margin, breaks = pretty(G$margin, n = N_BREAKS), plot = FALSE)
# freq = FALSE, or equidistant breaks would be plotted as counts
plot(h, freq = FALSE, col = ifelse(h$mids > 0, COL_FOC, COL_OTH), border = "white",
     main = "", xlab = "", ylab = "", axes = FALSE)
ax(1); ax(2, las = 1)
segments(0, 0, 0, max(h$density), lty = 2, col = "grey25", lwd = 3)
mtext(sprintf("%s z  -  best competing clade z", FOCAL_LABEL), side = 1, line = LAB_LINE,
      cex = LAB_CEX)
mtext(sprintf("density (n = %d)", nrow(G)), side = 2, line = YLAB_LINE, cex = LAB_CEX)
invisible(dev.off())

# summary of z and margin across maps (the % ahead is a frequency, not a p-value)
cat(sprintf("%s written (T_factor %.2f, %d maps)\n", OUT, T_USE, nrow(G)))
cat(sprintf("%s z: median %.2f | best competing clade z: median %.2f\n",
            FOCAL_LABEL, median(G$fz), median(G$best)))
cat(sprintf("margin: median %+.2f (95%% %+.2f to %+.2f); %s ahead in %.1f%% of maps\n",
            median(G$margin), quantile(G$margin, .025), quantile(G$margin, .975),
            FOCAL_LABEL, 100 * mean(G$margin > 0)))

# p_maxT from 15_timeslice.R (corrected for choosing the focal clade): spread across maps
if (file.exists(PERMAP)) {
  P <- fread(PERMAP)[T_factor == T_USE]
  cat(sprintf("p_maxT : median %.4f (95%% across maps %.4f to %.4f); %.1f%% of maps < 0.05\n",
              median(P$p_maxT), quantile(P$p_maxT, .025), quantile(P$p_maxT, .975),
              100 * mean(P$p_maxT < 0.05)))
  cat(sprintf("p_marg : median %.4f   (uncorrected, focal clade only)\n", median(P$p_marg)))
} else {
  cat(sprintf("p_maxT : %s not found -- rerun 15_timeslice.R to write it\n", PERMAP))
}
