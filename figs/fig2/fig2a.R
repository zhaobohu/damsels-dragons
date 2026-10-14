#!/usr/bin/env Rscript
# Fig. 2A: dated tree of every karyotyped species, branches coloured by haploid chromosome number

suppressPackageStartupMessages({
  library(ape); library(data.table); library(phytools); library(viridisLite)
})

TREE  <- "../../fusion_rates/tree/odon.dated.tre"
TAXA  <- "../../fusion_rates/taxa.tsv"
OUT   <- "a.bigtree.svg"
SHOW_LABELS  <- TRUE     # 426 species names; FALSE if you will annotate later
BAR_X        <- 0.02     # colour key left edge, fraction of plot width
BAR_Y        <- 0.06     # gap between tip labels and key, fraction of plot height
BAR_LEN      <- 0.16
BAR_CEX      <- 1.6
KEY_GAP      <- 0.03     # gap from the colour bar to the circle key
NEO_PCH       <- 21      # circle; bg is the fill, col the border
NEO_FILL      <- "#B22222"   # fill for species with a recorded sex-autosome fusion
NEO_OPEN      <- "#386cb0"     # fill for species without
NEO_BORDER    <- "grey25"
NEO_CEX       <- 0.9     # 426 tips are ~5 pt apart
NEO_LWD       <- 0.6
NEO_DY        <- 0.012   # circle drawn this far below the tip, fraction of plot height
WIDTH  <- 32     # inches
HEIGHT <- 7     # inches
CEX   <- 0.62    # tip labels
EDGE_LWD  <- 2      # branch thickness; near 10 the 426 branches merge
AXIS_LWD  <- 5      # axis line and ticks
AXIS_CEX  <- 2.4    # axis number size
AXIS_TCK  <- -0.035 # tick length, fraction of plot width (more negative = longer)
AXIS_NUM_LINE <- 3.0  # distance of the tick numbers from the axis; must exceed the tick length
AXIS_LINE <- 6.5    # distance of "Ma ago" from the axis
AXIS_TITLE_CEX <- 2.4

# read the tree and taxa table, keep the karyotyped species
tr <- read.tree(TREE); tr$tip.label <- sub("^_R_", "", tr$tip.label)
tx <- fread(TAXA); tx[, tip := gsub(" ", "_", species)]
tx <- tx[tip %chin% tr$tip.label]

# prune to every karyotyped species and ladderize
sk <- keep.tip(tr, tx$tip)
sk <- ladderize(sk, right = FALSE)
info <- tx[match(sk$tip.label, tip)]
sk$tip.label <- sub("_", " ", sk$tip.label)

nchr <- setNames(as.numeric(info$nchr), sk$tip.label)

# warn if the branches will be thicker than the gap between tips
tip_gap_pt <- (WIDTH - 2) * 72 / Ntip(sk)
if (EDGE_LWD > tip_gap_pt)
  warning(sprintf("EDGE_LWD %.1f exceeds the %.1f pt spacing between %d tips -- branches will merge",
                  EDGE_LWD, tip_gap_pt, Ntip(sk)), call. = FALSE)

# bottom margin from the measured width of the rotated tip labels
lab_ext_in <- if (SHOW_LABELS)
  max(strwidth(sk$tip.label, cex = CEX, font = 3, units = "inches")) else 0.15
BOT_MAR <- lab_ext_in / par("csi") + 1.2 + 4.5

# ancestral haploid chromosome number along the branches (takes a minute or two)
cat("fitting ancestral states on", Ntip(sk), "tips (this takes a minute)...\n")
cm <- setMap(contMap(sk, nchr, plot = FALSE), viridis(600))

# draw the tree downward; plot.contMap ignores par(mar), and legend = FALSE avoids a blank canvas
svg(OUT, width = WIDTH, height = HEIGHT)
plot(cm, direction = "downwards", lwd = EDGE_LWD, outline = FALSE, ftype = "off",
     legend = FALSE, mar = c(BOT_MAR, AXIS_LINE + 3.0, 0.6, 0.6))
par(xpd = NA)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xx <- lp$xx[seq_len(Ntip(sk))]; yy <- lp$yy[seq_len(Ntip(sk))]
sc <- diff(range(lp$yy[seq_len(Ntip(sk) + sk$Nnode)]))

# circle below every tip: filled where a fusion is recorded, open where not
neo <- info$sex_system != "X"
uu <- par("usr")
points(xx, yy - NEO_DY * (uu[4] - uu[3]), pch = NEO_PCH, cex = NEO_CEX,
       lwd = NEO_LWD, bg = ifelse(neo, NEO_FILL, NEO_OPEN), col = NEO_BORDER)

# tip labels by hand, rotated 90 (plot.phylo ignores srt)
if (SHOW_LABELS)
  text(xx, yy - NEO_DY * 2.2 * sc, sk$tip.label,
       srt = 90, adj = c(1, 0.5), cex = CEX, font = 3)

# time axis from the plotted coordinates (tips at y = 0, so y is age in Ma)
root_y <- lp$yy[Ntip(sk) + 1]
ax <- pretty(c(0, root_y), n = 4); ax <- ax[ax <= root_y]
axis(2, at = ax, labels = ax, cex.axis = AXIS_CEX, las = 1, col = "grey40",
     col.axis = "grey30", lwd = AXIS_LWD, tck = AXIS_TCK,
     mgp = c(3, AXIS_NUM_LINE, 0))
mtext("Ma ago", side = 2, line = AXIS_LINE, cex = AXIS_TITLE_CEX, col = "grey30")

# colour bar by hand, below the rotated labels (label width measured in inches, converted to y)
u <- par("usr")
y_per_in <- (u[4] - u[3]) / par("pin")[2]
lab_ext_u <- if (SHOW_LABELS)
  max(strwidth(sk$tip.label, cex = CEX, font = 3, units = "inches")) * y_per_in else 0
ky <- -lab_ext_u - BAR_Y * (u[4] - u[3])
kh <- 0.018 * (u[4] - u[3])
kx <- u[1] + BAR_X * (u[2] - u[1])
kw <- BAR_LEN * (u[2] - u[1])
px <- seq(kx, kx + kw, length.out = 121)
rect(px[-121], ky, px[-1], ky + kh, col = viridis(120), border = NA)
text(kx,      ky - 0.006 * (u[4] - u[3]), cm$lims[1], adj = c(0.5, 1), cex = BAR_CEX)
text(kx + kw, ky - 0.006 * (u[4] - u[3]), cm$lims[2], adj = c(0.5, 1), cex = BAR_CEX)
text(kx + kw / 2, ky + kh + 0.006 * (u[4] - u[3]), "haploid chromosomes",
     adj = c(0.5, 0), cex = BAR_CEX)

# circle key, to the right of the colour bar
w <- u[2] - u[1]
points(kx + kw + KEY_GAP * w, ky + kh / 2, pch = NEO_PCH, cex = NEO_CEX * 1.6,
       lwd = NEO_LWD, bg = NEO_FILL, col = NEO_BORDER)
text(kx + kw + (KEY_GAP + 0.012) * w, ky + kh / 2, "fusion recorded",
     adj = c(0, 0.5), cex = BAR_CEX * 0.9)
points(kx + kw + (KEY_GAP + 0.10) * w, ky + kh / 2, pch = NEO_PCH, cex = NEO_CEX * 1.6,
       lwd = NEO_LWD, bg = NEO_OPEN, col = NEO_BORDER)
text(kx + kw + (KEY_GAP + 0.112) * w, ky + kh / 2, "no fusion",
     adj = c(0, 0.5), cex = BAR_CEX * 0.9)
invisible(dev.off())
cat(sprintf("%s: %d species (%d on the tree); tip spacing %.1f pt, EDGE_LWD %.1f%s\n",
            OUT, Ntip(sk), Ntip(tr), tip_gap_pt, EDGE_LWD,
            sprintf("; haploid chromosomes %d-%d", min(nchr), max(nchr))))

# tip order left to right as plotted (by plotted x, not sk$tip.label order)
tip_order <- data.frame(
  position    = seq_len(Ntip(sk)),
  species     = sk$tip.label[order(xx)],
  family      = info$family[order(xx)],
  nchr        = info$nchr[order(xx)],
  sex_system  = info$sex_system[order(xx)],
  stringsAsFactors = FALSE)
write.csv(tip_order, sub("\\.svg$", ".tiporder.csv", OUT), row.names = FALSE)
cat(sprintf("tip order written to %s\n", sub("\\.svg$", ".tiporder.csv", OUT)))
print(tip_order, row.names = FALSE, max = 5 * nrow(tip_order))
