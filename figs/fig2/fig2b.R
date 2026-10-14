#!/usr/bin/env Rscript
# Fig. 2B: the Aeshna clade, branches coloured by reconstructed haploid chromosome number

suppressPackageStartupMessages({
  library(ape); library(data.table); library(phytools); library(viridisLite)
})

TREE   <- "../../fusion_rates/tree/odon.dated.tre"
TAXA <- "../../fusion_rates/taxa.tsv"
OUT  <- "b.aeshna.tree.svg"
GENUS  <- "Aeshna"

WIDTH   <- 32     # inches
HEIGHT  <- 8      # inches
LWD     <- 12     # branch thickness; scale with the canvas (~0.35 per inch)
FSIZE   <- 2      # tip label size
OFFSET  <- 5.5    # gap between tip and label, in character widths; must clear NEO_DY
NEO_PCH    <- 21        # circle; bg is the fill, col the border
NEO_FILL   <- "#B22222"   # species with a recorded sex-autosome fusion
NEO_OPEN   <- "#386cb0"   # species without
NEO_BORDER <- "grey25"
NEO_LWD    <- 1.5
NEO_CEX <- 2.5
NEO_DY  <- 0.035  # circle drawn this far below the tip (clear of the thick branch end),
                  # fraction of plot height
BAR_X   <- 0.62   # bar left edge, fraction of plot width
BAR_Y   <- 0.02   # bar height above the bottom, fraction of plot height
BAR_LEN <- 0.22   # bar length, fraction of plot width
BAR_CEX <- 1.6
KEY_X   <- 0.02   # circle key, left edge as a fraction of plot width
KEY_GAP <- 0.16   # gap to the second key entry
AXIS_LWD  <- 5    # time axis, as in fig2a.R
AXIS_CEX  <- 2.4
AXIS_TCK  <- -0.035
AXIS_NUM_LINE <- 3.0  # tick numbers; must exceed the tick length
AXIS_LINE <- 6.5
AXIS_TITLE_CEX <- 2.4
AXIS_MAR  <- 11.5   # left margin (lines) reserved for the axis

# read the tree and taxa table, extract the GENUS clade
tr <- read.tree(TREE); tr$tip.label <- sub("^_R_", "", tr$tip.label)
tx <- fread(TAXA); tx[, tip := gsub(" ", "_", species)]
keep <- intersect(tx[genus == GENUS, tip], tr$tip.label)
if (length(keep) < 3) stop("fewer than 3 ", GENUS, " tips on the tree", call. = FALSE)

cl <- keep.tip(extract.clade(tr, getMRCA(tr, keep)), keep)
cl <- ladderize(cl, right = FALSE)
info <- tx[match(cl$tip.label, tip)]
nchr <- setNames(as.numeric(info$nchr), cl$tip.label)
neo  <- info$sex_system != "X"

# scale the tree depth to 1 so LWD and FSIZE do not depend on clade age (age = y * DEPTH_MA)
DEPTH_MA <- max(node.depth.edgelength(cl))
cl$edge.length <- cl$edge.length / DEPTH_MA

# ancestral haploid chromosome number along the branches
cm <- contMap(cl, nchr, plot = FALSE)
cm <- setMap(cm, viridis(600))

# draw the tree downward; plot.contMap ignores par(mar), and legend = FALSE avoids a blank canvas
svg(OUT, width = WIDTH, height = HEIGHT)
plot(cm, direction = "downwards", lwd = LWD, outline = FALSE,
     ftype = "i", fsize = FSIZE, offset = OFFSET, legend = FALSE,
     mar = c(1.0, AXIS_MAR, 1.0, 1.0))
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
par(xpd = NA)

# time axis: ticks at round ages, converted back to the scaled y
ax_ma <- pretty(c(0, DEPTH_MA), n = 4); ax_ma <- ax_ma[ax_ma <= DEPTH_MA]
axis(2, at = ax_ma / DEPTH_MA, labels = ax_ma, cex.axis = AXIS_CEX, las = 1,
     col = "grey40", col.axis = "grey30", lwd = AXIS_LWD, tck = AXIS_TCK,
     mgp = c(3, AXIS_NUM_LINE, 0))
mtext("Ma ago", side = 2, line = AXIS_LINE, cex = AXIS_TITLE_CEX, col = "grey30")

# circle below every tip: filled where a fusion is recorded, open where not
u0 <- par("usr")
points(lp$xx[seq_len(Ntip(cl))],
       lp$yy[seq_len(Ntip(cl))] - NEO_DY * (u0[4] - u0[3]),
       pch = NEO_PCH, cex = NEO_CEX, lwd = NEO_LWD,
       bg = ifelse(neo, NEO_FILL, NEO_OPEN), col = NEO_BORDER)

# colour bar at bottom right
u <- par("usr")
add.color.bar(leg = BAR_LEN * (u[2] - u[1]), cols = cm$cols,
              title = "haploid chromosomes", lims = cm$lims, digits = 0,
              prompt = FALSE, subtitle = "", fsize = BAR_CEX, lwd = LWD,
              x = u[1] + BAR_X * (u[2] - u[1]),
              y = u[3] + BAR_Y * (u[4] - u[3]))

# circle key at bottom left
ky <- u[3] + BAR_Y * (u[4] - u[3]); w <- u[2] - u[1]
points(u[1] + KEY_X * w, ky, pch = NEO_PCH, cex = NEO_CEX,
       lwd = NEO_LWD, bg = NEO_FILL, col = NEO_BORDER)
text(u[1] + (KEY_X + 0.015) * w, ky, "fusion recorded",
     adj = c(0, 0.5), cex = BAR_CEX * 0.8)
points(u[1] + (KEY_X + KEY_GAP) * w, ky, pch = NEO_PCH, cex = NEO_CEX,
       lwd = NEO_LWD, bg = NEO_OPEN, col = NEO_BORDER)
text(u[1] + (KEY_X + KEY_GAP + 0.015) * w, ky, "no fusion",
     adj = c(0, 0.5), cex = BAR_CEX * 0.8)
invisible(dev.off())
cat(sprintf("%s: %d tips, %d neo-sex, haploid chromosomes %d-%d\n",
            OUT, Ntip(cl), sum(neo), min(nchr), max(nchr)))
