# Fig. 1A: riparian plots of the 29-genome GENESPACE run, one per order
# (needs genespace/state.rds from genespace.R; writes fig1a_odonata.svg, fig1a_ephemeroptera.svg)

library(GENESPACE)
library(data.table)
library(ggplot2)
library(svglite)

gs <- readRDS(file.path("genespace", "state.rds"))
gs$params$nCores <- 1L

# the 29 genomes of the run; the two panels below must cover them all
plotOrder <- c(
  # Anisoptera
  "Aeshna_juncea", "Anax_junius", "Anax_strenuus",          # Aeshnidae
  "Tanypteryx_hageni",                                      # Petaluridae
  "Cordulegaster_boltonii",                                 # Cordulegastridae
  "Somatochlora_meridionalis", "Somatochlora_arctica",       # Corduliidae
  "Orthetrum_coerulescens", "Sympetrum_striolatum",
  "Sympetrum_sanguineum", "Pantala_flavescens",             # Libellulidae
  # Zygoptera
  "Archilestes_grandis",                                    # Lestidae
  "Hetaerina_titia",                                        # Calopterygidae
  "Heteragrion_aequatoriale",                               # Heteragrionidae
  "Platycnemis_pennipes",                                   # Platycnemididae
  "Pyrrhosoma_nymphula", "Ceriagrion_tenellum",
  "Ischnura_elegans", "Ischnura_senegalensis",              # Coenagrionidae
  # Ephemeroptera
  "Serratella_ignita",                                      # Ephemerellidae
  "Ephemera_vulgata", "Ephemera_orientalis",                # Ephemeridae
  "Ecdyonurus_torrentis", "Parafronurus_youi",              # Heptageniidae
  "Siphlonurus_alternatus",                                 # Siphlonuridae
  "Baetis_rhodani", "Baetis_atlanticus",
  "Cloeon_dipterum", "Procloeon_bifidum"                    # Baetidae
)

# scaffold -> chromosome label (keys are the gff seqids); "Y" and "?" (unplaced) are not drawn
map <- c(
  # Aeshna juncea
  "OZ369542.1"="1","OZ369543.1"="2","OZ369544.1"="3","OZ369545.1"="4",
  "OZ369546.1"="5","OZ369547.1"="6","OZ369548.1"="7","OZ369549.1"="8",
  "OZ369550.1"="9","OZ369551.1"="10","OZ369552.1"="11","OZ369553.1"="12",
  "OZ369554.1"="X","OZ369555.1"="Y",
  # Anax junius
  "OZ349800.1"="1","OZ349801.1"="2","OZ349802.1"="3","OZ349803.1"="4",
  "OZ349804.1"="5","OZ349805.1"="6","OZ349806.1"="7","OZ349807.1"="8",
  "OZ349808.1"="9","OZ349809.1"="10","OZ349810.1"="11","OZ349811.1"="12",
  "OZ349812.1"="13","OZ349813.1"="X",
  # Anax strenuus
  "OZ282522.1"="1","OZ282523.1"="2","OZ282524.1"="3","OZ282525.1"="4",
  "OZ282526.1"="5","OZ282527.1"="6","OZ282528.1"="7","OZ282529.1"="8",
  "OZ282530.1"="9","OZ282531.1"="10","OZ282532.1"="11","OZ282533.1"="12",
  "OZ282534.1"="13","OZ282535.1"="X",
  # Archilestes grandis
  "OZ349685.1"="1","OZ349686.1"="2","OZ349687.1"="3","OZ349688.1"="4",
  "OZ349689.1"="5","OZ349690.1"="6","OZ349691.1"="7","OZ349692.1"="8",
  "OZ349693.1"="9","OZ349694.1"="10","OZ349695.1"="11","OZ349696.1"="X",
  # Baetis atlanticus
  "OZ373111.1"="1","OZ373112.1"="2","OZ373113.1"="3","OZ373114.1"="4",
  "OZ373115.1"="5",
  # Baetis rhodani
  "OZ391881.1"="1","OZ391882.1"="2","OZ391883.1"="3","OZ391884.1"="4",
  "OZ391885.1"="5",
  # Ceriagrion tenellum
  "OY720632.1"="1","OY720633.1"="2","OY720634.1"="3","OY720635.1"="4",
  "OY720636.1"="5","OY720637.1"="6","OY720639.1"="7","OY720640.1"="8",
  "OY720641.1"="9","OY720642.1"="10","OY720643.1"="11","OY720644.1"="12",
  "OY720645.1"="13","OY720638.1"="X",
  # Cloeon dipterum
  "OX451259.1"="1","OX451260.1"="2","OX451261.1"="3","OX451262.1"="4",
  "OX451263.1"="X",
  # Cordulegaster boltonii
  "OZ297306.1"="1","OZ297307.1"="2","OZ297308.1"="3","OZ297309.1"="4",
  "OZ297310.1"="5","OZ297311.1"="6","OZ297312.1"="7","OZ297313.1"="8",
  "OZ297314.1"="9","OZ297315.1"="10","OZ297316.1"="11","OZ297317.1"="12",
  "OZ297318.1"="X",
  # Ecdyonurus torrentis
  "OX439128.1"="1","OX439129.1"="2","OX439130.1"="3","OX439131.1"="4",
  "OX439132.1"="5","OX439133.1"="6","OX439134.1"="7","OX439135.1"="8",
  "OX439136.1"="9","OX439137.1"="10","OX439127.1"="X",
  # Ephemera orientalis
  "CM180864.1"="1","CM180865.1"="2","CM180866.1"="3","CM180867.1"="4",
  "CM180868.1"="5","CM180869.1"="X",
  # Ephemera vulgata
  "OZ367565.1"="1","OZ367566.1"="2","OZ367567.1"="3","OZ367568.1"="4",
  "OZ367569.1"="5","OZ367570.1"="6",
  # Hetaerina titia
  "CM073768.1"="1","CM073769.1"="2","CM073770.1"="3","CM073771.1"="4",
  "CM073772.1"="5","CM073773.1"="6","CM073774.1"="7","CM073775.1"="8",
  "CM073776.1"="9","CM073777.1"="10","CM073778.1"="11","CM073779.1"="X",
  "JAVFHJ010000013.1"="?","JAVFHJ010000014.1"="?",
  # Heteragrion aequatoriale
  "OZ485862.1"="1","OZ485863.1"="2","OZ485864.1"="3","OZ485865.1"="4",
  "OZ485866.1"="5","OZ485867.1"="6","OZ485868.1"="7","OZ485869.1"="8",
  "OZ485870.1"="9","OZ485871.1"="10","OZ485872.1"="11","OZ485873.1"="12",
  "OZ485874.1"="X",
  # Ischnura elegans
  "NC_060246.1"="1","NC_060247.1"="2","NC_060248.1"="3","NC_060249.1"="4",
  "NC_060250.1"="5","NC_060251.1"="6","NC_060252.1"="7","NC_060253.1"="8",
  "NC_060254.1"="9","NC_060255.1"="10","NC_060256.1"="11","NC_060257.1"="12",
  "NC_060258.1"="13","NC_060259.1"="X","NW_025791657.1"="?","NW_025791658.1"="?",
  "NW_025791659.1"="?","NW_028141064.1"="?",
  # Ischnura senegalensis
  "OZ176221.1"="1","OZ176222.1"="2","OZ176223.1"="3","OZ176224.1"="4",
  "OZ176225.1"="5","OZ176226.1"="6","OZ176227.1"="7","OZ176228.1"="8",
  "OZ176229.1"="9","OZ176231.1"="10","OZ176232.1"="11","OZ176233.1"="12",
  "OZ176234.1"="13","OZ176230.1"="X",
  # Orthetrum coerulescens
  "OZ294519.1"="1","OZ294520.1"="2","OZ294521.1"="3","OZ294522.1"="4",
  "OZ294523.1"="5","OZ294524.1"="6","OZ294525.1"="7","OZ294526.1"="8",
  "OZ294527.1"="9","OZ294528.1"="10","OZ294529.1"="11","OZ294530.1"="12",
  "OZ294531.1"="X","OZ294532.1"="Y",
  # Pantala flavescens
  "CM036878.1"="1","CM036879.1"="2","CM036880.1"="3","CM036881.1"="4",
  "CM036882.1"="5","CM036883.1"="6","CM036884.1"="7","CM036885.1"="8",
  "CM036886.1"="9","CM036887.1"="10","CM036888.1"="11","CM036889.1"="X",
  # Parafronurus youi
  "CM104470.1"="1","CM104471.1"="2","CM104472.1"="3","CM104473.1"="4",
  "CM104474.1"="5","CM104475.1"="6","CM104476.1"="7","CM104477.1"="8",
  "CM104478.1"="9","CM104479.1"="10","CM104480.1"="11",
  # Platycnemis pennipes
  "OW121847.1"="1","OW121848.1"="2","OW121849.1"="3","OW121850.1"="4",
  "OW121851.1"="5","OW121852.1"="6","OW121853.1"="7","OW121854.1"="8",
  "OW121855.1"="9","OW121856.1"="10","OW121857.1"="11","OW121858.1"="12",
  "OW121846.1"="X",
  # Procloeon bifidum
  "OZ259029.1"="1","OZ259030.1"="2","OZ259031.1"="3","OZ259032.1"="4",
  "OZ259033.1"="5",
  # Pyrrhosoma nymphula
  "OY754327.1"="1","OY754328.1"="2","OY754329.1"="3","OY754330.1"="4",
  "OY754331.1"="5","OY754333.1"="6","OY754334.1"="7","OY754335.1"="8",
  "OY754336.1"="9","OY754337.1"="10","OY754338.1"="11","OY754339.1"="12",
  "OY754340.1"="13","OY754332.1"="X",
  # Serratella ignita
  "OZ482020.1"="1","OZ482021.1"="2","OZ482022.1"="3","OZ482023.1"="4",
  "OZ482024.1"="5","OZ482025.1"="6","OZ482026.1"="7","OZ482027.1"="X",
  # Siphlonurus alternatus
  "OX463778.1"="1","OX463779.1"="2","OX463780.1"="3","OX463781.1"="4",
  "OX463782.1"="5","OX463783.1"="6","OX463784.1"="7","OX463785.1"="8",
  "OX463786.1"="9","OX463787.1"="10","OX463788.1"="X",
  # Somatochlora arctica
  "OZ407513.1"="1","OZ407514.1"="2","OZ407515.1"="3","OZ407516.1"="4",
  "OZ407517.1"="5","OZ407518.1"="6","OZ407519.1"="7","OZ407520.1"="8",
  "OZ407521.1"="9","OZ407522.1"="10","OZ407523.1"="11","OZ407524.1"="12",
  "OZ407525.1"="X",
  # Somatochlora meridionalis
  "OZ317805.1"="1","OZ317806.1"="2","OZ317807.1"="3","OZ317808.1"="4",
  "OZ317809.1"="5","OZ317810.1"="6","OZ317811.1"="7","OZ317812.1"="8",
  "OZ317813.1"="9","OZ317814.1"="10","OZ317815.1"="11","OZ317816.1"="12",
  # Sympetrum sanguineum
  "OZ078875.1"="1","OZ078876.1"="2","OZ078877.1"="3","OZ078878.1"="4",
  "OZ078879.1"="5","OZ078880.1"="6","OZ078881.1"="7","OZ078883.1"="8",
  "OZ078884.1"="9","OZ078885.1"="10","OZ078886.1"="11","OZ078887.1"="12",
  "OZ078882.1"="X",
  # Sympetrum striolatum
  "OX388341.1"="1","OX388342.1"="2","OX388343.1"="3","OX388344.1"="4",
  "OX388345.1"="5","OX388346.1"="6","OX388347.1"="7","OX388348.1"="8",
  "OX388349.1"="9","OX388350.1"="10","OX388351.1"="11","OX388352.1"="X",
  # Tanypteryx hageni
  "CM054112.1"="1","CM054113.1"="2","CM054114.1"="3","CM054115.1"="4",
  "CM054116.1"="5","CM054117.1"="6","CM054118.1"="7","CM054119.1"="8",
  "CM054120.1"="9"
)

chrLabFun <- function(x) ifelse(x %in% names(map), unname(map[x]), x)
drop_chr  <- names(map)[map %in% c("?", "Y")]   # hide unplaced + Y scaffolds

pal    <- colorRampPalette(c("darkorange","skyblue","darkblue",
                             "purple","darkred","salmon"))
addThm <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))

# phased blocks for a reference: read <ref>_phasedBlks.csv, or compute it with
# plot_riparian() and save it (GENESPACE 1.3.1 plot_riparian does not write it)
get_phased <- function(ref, ids) {
  f <- file.path(gs$paths$riparian, sprintf("%s_phasedBlks.csv", ref))
  if (file.exists(f)) return(data.table::fread(f))
  rip <- GENESPACE::plot_riparian(
    gsParam = gs, refGenome = ref, genomeIDs = ids, useRegions = TRUE,
    pdfFile = tempfile(fileext = ".pdf"))
  data.table::fwrite(rip$blks, f)
  data.table::data.table(rip$blks)
}

# chromosomes to draw reversed: starting from the reference row (never flipped),
# flip a chromosome when its anchors, signed by block orientation and by the
# partner's flip, mostly run backwards relative to the already-settled neighbour row
auto_invert <- function(blk, order, ref = order[1], verbose = TRUE) {
  blk  <- data.table::as.data.table(blk)
  flip <- data.table::data.table(genome = character(), chr = character())
  k <- match(ref, order)
  if (is.na(k)) stop("auto_invert: reference ", ref, " is not in the row order")
  n <- length(order)
  steps <- c(if (k < n) (k + 1):n, if (k > 1) (k - 1):1)
  for (i in steps) {
    up <- order[if (i > k) i - 1 else i + 1]   # neighbour on the reference side
    g  <- order[i]
    b <- blk[genome1 == up & genome2 == g,
             .(chrUp = chr1, chr = chr2, orient, w = pmin(nHits1, nHits2))]
    if (!nrow(b))                              # blocks stored the other way round
      b <- blk[genome1 == g & genome2 == up,
               .(chrUp = chr2, chr = chr1, orient, w = pmin(nHits1, nHits2))]
    if (!nrow(b)) next
    b[, sUp := ifelse(paste(up, chrUp) %in% paste(flip$genome, flip$chr), -1, 1)]
    v <- b[, .(score = sum(ifelse(orient == "+", 1, -1) * sUp * w)), by = chr]
    if (any(v$score < 0))
      flip <- rbind(flip, data.table::data.table(genome = g, chr = v[score < 0, chr]))
  }
  if (verbose && nrow(flip)) {
    cat("[flip] chromosomes drawn reversed:\n")
    print(flip[, .(chrs = paste(chrLabFun(chr), collapse = " ")), by = genome])
  }
  if (nrow(flip)) as.data.frame(flip) else NULL
}

# remove the "*" riparian_engine() appends to flipped chromosome labels
unstar <- function(pl) {
  for (i in seq_along(pl$layers)) {
    d <- pl$layers[[i]]$data
    if (is.data.frame(d) && "chrLab" %in% names(d)) {
      d <- data.table::copy(d)
      d$chrLab <- sub("\\*$", "", d$chrLab)
      pl$layers[[i]]$data <- d
    }
  }
  pl
}

# gene positions of every genome, without the unplaced and Y scaffolds
bed  <- GENESPACE::read_combBed(file.path(gs$paths$results, "combBed.txt"))
bed  <- bed[!(chr  %in% drop_chr)]

CHR_LAB_SIZE <- 10
CHR_EXPAND   <- 40 / CHR_LAB_SIZE - 2
stopifnot(CHR_EXPAND >= 0)

# panels: each is phased to its own reference, so braid colours are
# A. juncea chromosomes in the odonate panel and Ecdyonurus chromosomes in the mayfly panel
EPH <- c("Serratella_ignita", "Ephemera_vulgata", "Ephemera_orientalis",
         "Ecdyonurus_torrentis", "Parafronurus_youi", "Siphlonurus_alternatus",
         "Baetis_rhodani", "Baetis_atlanticus", "Cloeon_dipterum",
         "Procloeon_bifidum")

# rows top to bottom, as in the species tree (rooted between the orders, ladderized)
odoOrder <- c(
  "Archilestes_grandis", "Heteragrion_aequatoriale", "Hetaerina_titia",
  "Platycnemis_pennipes", "Ischnura_senegalensis", "Ischnura_elegans",
  "Pyrrhosoma_nymphula", "Ceriagrion_tenellum",                 # Zygoptera
  "Aeshna_juncea", "Anax_strenuus", "Anax_junius",              # Aeshnidae (reference)
  "Tanypteryx_hageni", "Cordulegaster_boltonii",
  "Somatochlora_meridionalis", "Somatochlora_arctica",
  "Sympetrum_striolatum", "Sympetrum_sanguineum",
  "Pantala_flavescens", "Orthetrum_coerulescens")
ephOrder <- c(
  "Procloeon_bifidum", "Cloeon_dipterum",
  "Baetis_rhodani", "Baetis_atlanticus",                        # Baetidae
  "Serratella_ignita", "Ephemera_vulgata", "Ephemera_orientalis",
  "Siphlonurus_alternatus", "Parafronurus_youi",
  "Ecdyonurus_torrentis")                                       # reference
stopifnot(setequal(ephOrder, EPH), setequal(odoOrder, setdiff(plotOrder, EPH)),
          anyDuplicated(c(odoOrder, ephOrder)) == 0)

PANELS <- list(
  list(ref = "Aeshna_juncea",        order = odoOrder, file = "fig1a_odonata.svg"),
  list(ref = "Ecdyonurus_torrentis", order = ephOrder, file = "fig1a_ephemeroptera.svg"))

# shared geometry of the two panels
IN_PER_ROW  <- 1.1  # distance between genome rows, inches
PANEL_W     <- 15   # plotting area width, inches (same x range in both panels)

# pass 1: phase, flip and build each panel on a scratch device of the final
# height, so chromosome bars get the same thickness in both
built <- lapply(PANELS, function(p) {
  ids  <- rev(p$order)                        # drawn bottom-up, so p$order reads top to bottom
  blkP <- get_phased(p$ref, ids)
  blkP <- blkP[!(chr1 %in% drop_chr | chr2 %in% drop_chr)]
  invP <- auto_invert(blkP, p$order, ref = p$ref)   # per-panel flips
  grDevices::pdf(NULL, width = PANEL_W, height = (length(ids) + 0.5) * IN_PER_ROW)
  pd <- GENESPACE::riparian_engine(
    blk               = blkP,
    bed               = bed,
    refGenome         = p$ref,
    genomeIDs         = ids,
    labelTheseGenomes = ids,
    palette           = pal,
    chrLabFun         = chrLabFun,
    chrBorderLwd      = 0.1,
    chrExpand         = CHR_EXPAND,
    chrFill           = "white",
    howSquare         = 5,
    chrLabFontSize    = CHR_LAB_SIZE,
    chrBorderCol      = "black",
    syntenyWeight     = 1,
    invertTheseChrs   = invP,
    addThemes         = addThm
  )
  grDevices::dev.off()
  c(p, list(n = length(ids), plot = pd$ggplotObj,
            chrs = data.table::as.data.table(pd$sourceData$chromosomes)))
})

# shared x range (genomes are centred on 0) and scale bar length
# (30% of the smallest genome, one significant digit, over both panels)
xhalf <- 1.005 * max(sapply(built, function(b) max(abs(c(b$chrs$x1, b$chrs$x2)))))
SCALE_GENES <- signif(0.3 * min(sapply(built, function(b)
    min(tapply(b$chrs$length, b$chrs$genome, sum)))), 1)
stopifnot(SCALE_GENES < xhalf)

# pass 2: replace GENESPACE's scale bar with the shared one, fix the axis ranges
# and set the plotting area in inches (svglite scratch device, so text is measured
# with the final fonts)
svglite::svglite(tempfile(fileext = ".svg"))
grobs <- lapply(built, function(b) {
  cph  <- max(b$chrs$y1 - b$chrs$y2)          # chromosome bar height, in rows
  y0   <- b$n + 1.5 * cph                     # bar sits above the top row, as in GENESPACE
  ylim <- c(1 - cph / 2 - 0.1, y0 + cph + 0.1)
  pl   <- unstar(b$plot)                 # no "*" on flipped chromosomes
  isBar <- vapply(pl$layers, function(l)
    is.data.frame(l$data) && "line" %in% names(l$data), logical(1))
  stopifnot(sum(isBar) == 2)                  # GENESPACE's bar: segments + label
  pl$layers <- pl$layers[!isBar]
  pl <- pl +
    ggplot2::annotate("segment",
      x = c(0, SCALE_GENES, 0), xend = c(0, SCALE_GENES, SCALE_GENES),
      y = c(y0, y0, y0 + cph / 2), yend = c(y0 + cph, y0 + cph, y0 + cph / 2),
      linewidth = 0.1, colour = "black") +
    ggplot2::annotate("text", x = SCALE_GENES / 2, y = y0 + cph / 2,
      label = paste(format(SCALE_GENES, big.mark = ",", scientific = FALSE), "genes"),
      size = CHR_LAB_SIZE * 0.36, vjust = 1.5) +
    ggplot2::coord_cartesian(xlim = c(-xhalf, xhalf), ylim = ylim, expand = FALSE)
  g   <- ggplot2::ggplotGrob(pl)
  pan <- g$layout[g$layout$name == "panel", ]
  g$widths[pan$l]  <- grid::unit(PANEL_W, "in")
  g$heights[pan$t] <- grid::unit(diff(ylim) * IN_PER_ROW, "in")
  g
})
wIn <- lapply(grobs, function(g) grid::convertWidth(g$widths, "in", valueOnly = TRUE))
hIn <- lapply(grobs, function(g) grid::convertHeight(g$heights, "in", valueOnly = TRUE))
grDevices::dev.off()
stopifnot(length(unique(lengths(wIn))) == 1)
wAll <- do.call(pmax, wIn)                    # same margins in both files

# write each panel with the shared margins
for (i in seq_along(built)) {
  g <- grobs[[i]]
  g$widths <- grid::unit(wAll, "in")
  svglite::svglite(built[[i]]$file, width = sum(wAll), height = sum(hIn[[i]]))
  grid::grid.newpage()
  grid::grid.draw(g)
  grDevices::dev.off()
  message(sprintf("[wrote] %s  (reference %s, %d genomes, %.2f x %.2f in, scale bar %s genes)",
                  built[[i]]$file, built[[i]]$ref, built[[i]]$n,
                  sum(wAll), sum(hIn[[i]]), format(SCALE_GENES, big.mark = ",")))
}
