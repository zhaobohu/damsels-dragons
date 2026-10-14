#!/usr/bin/env Rscript
# p_maxT per map, one beeswarm per slice time (Supplementary Results; not a manuscript figure)

suppressPackageStartupMessages({ library(data.table); library(devEMF) })

PERMAP <- "../../fusion_rates/timeslice_permap.tsv"
OUT    <- "cut_age.emf"
COLUMNS <- c(1, 1.25, 1.5, 2)        # T factors to draw, in order; others skipped

WIDTH <- 12; HEIGHT <- 9
SW   <- 0.42      # half-width of the swarm
STEP <- 0.013     # max sideways spacing between adjacent points in a row
NMAX <- 1200      # points drawn per column; larger sets are subsampled
W    <- 0.18      # box half-width
PT_CEX <- 0.6; PT_ALPHA <- 0.55
COLS <- c("#D85A30", "#E08A6A", "#B4B2A9", "#5F5E5A")   # one per column
YLIM_MIN <- NA    # y axis floor; NA = set from the data (smallest p is 1/(N_NULL + 1))
AXIS_LWD <- 3; AXIS_CEX <- 1.6; LAB_CEX <- 1.8; LAB_LINE <- 4.2; TOP_CEX <- 1.1

# read p_maxT, set the log10 y axis floor, clamp values below it
P <- fread(PERMAP)
P <- P[T_factor %in% COLUMNS]
floor_p <- if (is.na(YLIM_MIN)) 10^floor(log10(min(P$p_maxT))) else YLIM_MIN
YLIM <- c(floor_p, 1)
n_clamped <- sum(P$p_maxT < YLIM[1])
if (n_clamped) warning(sprintf("%d maps below the axis floor %.0e are drawn at the floor",
                               n_clamped, YLIM[1]),
                       call. = FALSE)
P[, y := log10(pmax(p_maxT, YLIM[1]))]

# x-axis labels: T factor and median slice age in Ma
med_age <- c()
for (tf in COLUMNS) med_age <- c(med_age, median(P[T_factor == tf, T_Ma]))
labels <- sprintf("%sx\n%.0f Ma", sub("\\.?0+$", "", sprintf("%.2f", COLUMNS)), med_age)

ylo <- log10(YLIM[1]); yhi <- log10(YLIM[2])

# empty plot with the p = 0.05 guide and both axes
set.seed(1)
emf(OUT, width = WIDTH, height = HEIGHT)
par(mar = c(LAB_LINE + 2.5, LAB_LINE + 3, 4.5, 1), xpd = NA)
plot(NA, xlim = c(0.4, length(COLUMNS) + 0.6), ylim = c(ylo, yhi), axes = FALSE,
     xlab = "", ylab = "")
abline(h = log10(0.05), lty = 3, col = "grey55", lwd = 2)
text(par("usr")[1], log10(0.05), "  p = 0.05", adj = c(0, -0.4), cex = 0.9, col = "grey45")
ticks <- c(10^seq(log10(YLIM[1]), -1), 0.05, 0.5, 1)
ticks <- sort(unique(ticks[ticks >= YLIM[1]]))
tick_lab <- formatC(ticks, format = "fg", digits = 4, drop0trailing = TRUE)  # not 1e-04
axis(2, at = log10(ticks), labels = tick_lab, las = 1, lwd = AXIS_LWD, cex.axis = AXIS_CEX,
     col = "grey25", col.axis = "grey15")
axis(1, at = seq_along(COLUMNS), labels = labels, tick = FALSE, line = 0.8,
     cex.axis = AXIS_CEX, padj = 0.5, col.axis = "grey15")
mtext("p_maxT per map", side = 2, line = LAB_LINE, cex = LAB_CEX)
mtext("slice time (multiple of focal crown age; Ma)", side = 1, line = LAB_LINE + 1.4,
      cex = LAB_CEX)

for (i in seq_along(COLUMNS)) {
  v <- P[T_factor == COLUMNS[i], y]

  # swarm: subsample, bin points into 90 rows, spread each row sideways
  pts <- v
  if (length(pts) > NMAX) pts <- sample(pts, NMAX)
  pts <- pts[is.finite(pts)]
  bin <- round((pts - ylo) / ((yhi - ylo) / 90))
  xoff <- rep(0, length(pts))
  for (b in unique(bin)) {
    k <- which(bin == b); k <- k[order(pts[k])]
    step <- min(STEP, 2 * SW / length(k))
    xoff[k] <- (seq_along(k) - (length(k) + 1) / 2) * step
  }
  points(i + xoff, pts, pch = 16, cex = PT_CEX, col = adjustcolor(COLS[i], PT_ALPHA))

  # box over all values: whiskers at 1.5 IQR, translucent IQR box, heavy median bar
  q <- quantile(v, c(.25, .5, .75)); iqr <- q[3] - q[1]
  bx <- c(min(v[v >= q[1] - 1.5 * iqr]), q[1], q[2], q[3], max(v[v <= q[3] + 1.5 * iqr]))
  segments(i, bx[1], i, bx[2]); segments(i, bx[4], i, bx[5])
  segments(i - W/2, bx[1], i + W/2, bx[1]); segments(i - W/2, bx[5], i + W/2, bx[5])
  rect(i - W, bx[2], i + W, bx[4], col = adjustcolor("white", 0.65), border = "black")
  segments(i - W, bx[3], i + W, bx[3], lwd = 2.6)

  # median p_maxT, % below 0.05 and n above the column
  pv <- P[T_factor == COLUMNS[i], p_maxT]
  top <- sprintf("median %.3f\n%.0f%% < 0.05", median(pv), 100 * mean(pv < 0.05))
  text(i, yhi + 0.06 * (yhi - ylo), top, cex = TOP_CEX, col = "grey15", font = 2,
       adj = c(0.5, 0))
  text(i, yhi + 0.01 * (yhi - ylo), sprintf("n = %d", length(v)), cex = 0.85, col = "grey25",
       adj = c(0.5, 0))
}
invisible(dev.off())

# summary table across cut ages
cat(sprintf("%s written; y axis floor %.0e (%d maps below it)\n", OUT, YLIM[1], n_clamped))
print(P[, .(maps = .N, median = round(median(p_maxT), 4), q025 = round(quantile(p_maxT, .025), 4),
            q975 = round(quantile(p_maxT, .975), 4), frac_sig = round(mean(p_maxT < 0.05), 3)),
        by = T_factor][order(T_factor)], row.names = FALSE)
