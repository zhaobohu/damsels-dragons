#!/usr/bin/env Rscript
# Fig. 2C: sex-autosome fraction among fusions per map: smoothed histograms, 95% interval, median

suppressPackageStartupMessages(library(data.table))

IN   <- "../../fusion_rates/fig_levels_permap.tsv"
OUT  <- "c.shapes.svg"
WIDTH  <- 32     # inches
AXIS_LWD   <- 5      # x axis line
AXIS_CEX   <- 2.4    # x axis numbers
AXIS_TCK   <- -0.035 # tick length (more negative = longer)
AXIS_NUM_LINE <- 3.0 # tick numbers; must exceed the tick length
LAB_CEX    <- 2.6    # x axis title
BAR_LWD    <- 13     # 95% interval bar
DOT_CEX    <- 3.6    # median dot
DOT_LWD    <- 4
CURVE_LWD  <- 4      # outline of each density
HEIGHT <- 12     # inches
FOCAL_LABEL <- "Aeshna"
BW <- 0.035      # one fixed bandwidth for all curves (the focal fraction takes few values)

# one entry per distribution, top to bottom, with its colour
P <- fread(IN)
series <- list(
  list(lab = "random expectation", x = P$null_frac,                     col = "#B4B2A9"),
  list(lab = FOCAL_LABEL,          x = P$foc_frac[!is.na(P$foc_frac)], col = "#B22222"),
  list(lab = "rest of the tree",   x = P$bg_frac[!is.na(P$bg_frac)],   col = "#386cb0"))

# density of each series, the tallest peak, and the shared x range
dens <- list(); ymax <- -Inf; xq <- c()
for (k in seq_along(series)) {
  dens[[k]] <- density(series[[k]]$x, bw = BW, from = 0, to = 1, n = 512)
  ymax <- max(ymax, dens[[k]]$y)
  xq <- c(xq, quantile(series[[k]]$x, c(.001, .999)))
}
xr   <- range(xq)
xr   <- c(max(0, xr[1] - 0.03), min(1, xr[2] + 0.03))

# empty plot with the x axis
svg(OUT, width = WIDTH, height = HEIGHT)
par(mar = c(9.5, 1.2, 1.0, 1.2), xpd = NA, mgp = c(2.2, 2.2, 0))
ROW <- 1; GAP <- 0.08            # row height and gap between rows
plot(NA, xlim = xr, ylim = c(0, 3 * (ROW + GAP) + 0.05),
     yaxt = "n", xaxt = "n", ylab = "", xlab = "", bty = "n")
axis(1, lwd = AXIS_LWD, cex.axis = AXIS_CEX, tck = AXIS_TCK, col = "grey25",
     col.axis = "grey15", mgp = c(3, AXIS_NUM_LINE, 0))
mtext("sex-autosome fusions / all fusions", side = 1, line = 6.0, cex = LAB_CEX)

# one row per series: density on top, 95% interval bar and median dot below
for (k in 3:1) {
  base <- (3 - k) * (ROW + GAP)          # Aeshna at the bottom, null at the top
  d <- dens[[k]]; x <- series[[k]]$x
  col <- series[[k]]$col
  yy <- base + 0.22 + 0.72 * ROW * d$y / ymax
  polygon(c(d$x, rev(d$x)), c(yy, rep(base + 0.22, length(d$x))),
          col = paste0(col, "99"), border = col, lwd = CURVE_LWD)
  q <- quantile(x, c(.025, .5, .975))
  segments(q[1], base + 0.1, q[3], base + 0.1, col = col, lwd = BAR_LWD, lend = 1)
  points(q[2], base + 0.1, pch = 21, bg = "white", col = col, cex = DOT_CEX, lwd = DOT_LWD)
}
invisible(dev.off())

# print median and 95% interval of each series
for (k in seq_along(series)) {
  q <- quantile(series[[k]]$x, c(.025, .5, .975))
  cat(sprintf("  %-20s median %.3f  (%.3f - %.3f)  n = %d\n",
              series[[k]]$lab, q[2], q[1], q[3], length(series[[k]]$x)))
}
cat("written:", OUT, "\n")
