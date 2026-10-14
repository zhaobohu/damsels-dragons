# Fig. 1B (% orthologs retained on the homologous chromosome) and 1C (ARI) against divergence
library(data.table)

# paths: absolute or relative to this script's directory; a command-line argument overrides PAIRS
DATA    <- "../../palaeoptera_synteny/pairwise_synteny"
OUTDIR  <- "."                        # "." = this script's own directory
PAIRS   <- "pairs_o2o.tsv"            # inside DATA
OUT_B   <- "fig1b_retained"           # .svg added when written
OUT_C   <- "fig1c_ari"

# appearance, one entry per panel
LEG_LAB     <- c("Odonata", "Ephemeroptera")
XLAB        <- "substitutions / site"  # patristic distance on the IQ-TREE species tree
PANELS <- list(
  list(file = OUT_B, y = "retained", ylim = c(0, 100), leg = "bottomleft",
       ylab = "% orthologs on homologous chromosome"),
  list(file = OUT_C, y = "ari",      ylim = c(0, 1),   leg = "bottomleft",
       ylab = "chromosome-assignment agreement (ARI)"))

WIDTH <- 15; HEIGHT <- 14             # inches
CEX_LAB <- 3; CEX_AXIS <- 3; CEX_LEG <- 3; CEX_PT <- 3; TCL <- -1.5
col_odo  <- adjustcolor("#1f78b4", 0.4)
col_eph  <- adjustcolor("#b2df8a", 0.4)

# resolve paths against the script's own directory
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
} else getwd()
resolve <- function(p) normalizePath(
  if (grepl("^(/|~)", p)) path.expand(p) else file.path(script_dir, p), mustWork = FALSE)

arg <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(arg)) PAIRS <- arg
PAIRS  <- if (file.exists(PAIRS)) PAIRS else file.path(resolve(DATA), PAIRS)
OUTDIR <- resolve(OUTDIR)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
G_ODO <- "Odo-Odo"; G_EPH <- "Eph-Eph"

# margins (lines) to fit tick labels and axis titles at these text sizes
MGP <- c(abs(TCL) + 1.9 * CEX_AXIS + 0.6, abs(TCL) + 0.7 * CEX_AXIS, 0)
MAR <- c(MGP[1] + 1.6 * CEX_LAB, MGP[1] + 1.6 * CEX_LAB, 2,
         1.1 * CEX_AXIS)

# read pairs, keep within-order pairs only
d <- fread(PAIRS)
miss <- setdiff(c("class", "dist", "ari", "nest"), names(d))
if (length(miss)) stop(basename(PAIRS), " is missing: ", paste(miss, collapse = ", "))
d[, retained  := 100 * nest]
d[, displaced := 100 - retained]      # console cross-check with 03_stats.R only
w <- d[class %chin% c(G_ODO, G_EPH)]
if (!nrow(w)) stop("no within-order pairs")

# medians within the shared divergence range, to cross-check with 03_stats.R
lim <- w[class == G_ODO, max(dist)]   # shared-range bound
s   <- w[dist <= lim]
cat(sprintf("[shared range] dist <= %.4f   %d of %d mayfly pairs beyond it\n",
            lim, w[class == G_EPH & dist > lim, .N], w[class == G_EPH, .N]))
cat("[check - must match 03_stats.R section 2]\n")
print(s[, .(pairs = .N, displaced = round(median(displaced), 2),
            retained = round(median(retained), 2), ari = round(median(ari), 3)),
        by = class][order(class)])

# x-axis from 0 to the first labelled tick at or beyond the deepest pair
XTICKS <- pretty(c(0, max(w$dist)))
XLIM   <- range(XTICKS)

# one svg per panel: points only, axes meeting at (0, 0)
for (p in PANELS) {
  below <- sum(w[[p$y]] < p$ylim[1])
  if (below) message(sprintf("[note] %s: %d pair(s) below %g (min %.3f), drawn under the x-axis",
                             p$y, below, p$ylim[1], min(w[[p$y]])))

  out <- paste0(file.path(OUTDIR, p$file), ".svg")
  svglite::svglite(out, width = WIDTH, height = HEIGHT, bg = "transparent")
  par(mar = MAR, cex.lab = CEX_LAB, cex.axis = CEX_AXIS, tcl = TCL, mgp = MGP)

  plot(NA, xlim = XLIM, ylim = p$ylim, xaxs = "i", yaxs = "i", axes = FALSE,
       xlab = XLAB, ylab = p$ylab)
  axis(1, at = XTICKS)
  axis(2, at = pretty(p$ylim), las = 1)

  # xpd: points on a limit (e.g. 100% retained) are drawn whole, not clipped
  op <- par(xpd = NA)
  e <- w[class == G_EPH]; o <- w[class == G_ODO]
  points(e$dist, e[[p$y]], pch = 16, cex = CEX_PT, col = col_eph)
  points(o$dist, o[[p$y]], pch = 16, cex = CEX_PT, col = col_odo)
  par(op)

  legend(p$leg, bty = "n", pch = 16, pt.cex = CEX_PT, cex = CEX_LEG,
         col = c(col_odo, col_eph), legend = LEG_LAB)
  dev.off()
  message("[wrote] ", out)
}
