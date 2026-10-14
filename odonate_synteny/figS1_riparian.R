# Fig. S1: riparian plot of the 10-odonate GENESPACE run (needs genespace/state.rds from genespace.R)

library(GENESPACE)
library(data.table)
library(ggplot2)
library(svglite)

gs <- readRDS(file.path("genespace", "state.rds"))
gs$params$nCores <- 1L

refGenome <- "Aeshna_juncea"

# stack order; genomes are passed to GENESPACE as rev(plotOrder)
plotOrder <- c(
  "Archilestes_grandis", "Platycnemis_pennipes", "Ischnura_elegans", "Cordulegaster_boltonii",
  "Somatochlora_arctica", "Sympetrum_striolatum", "Orthetrum_coerulescens",
  "Pantala_flavescens", "Anax_junius", "Aeshna_juncea"
)
genomeIDs <- rev(plotOrder)

# scaffold -> chromosome label (keys are the gff seqids); "Y" and "?" (unplaced) are not drawn
map <- c(
  # Aeshna juncea  (X = OZ369554.1, Y = OZ369555.1)
  "OZ369542.1"="1","OZ369543.1"="2","OZ369544.1"="3","OZ369545.1"="4",
  "OZ369546.1"="5","OZ369547.1"="6","OZ369548.1"="7","OZ369549.1"="8",
  "OZ369550.1"="9","OZ369551.1"="10","OZ369552.1"="11","OZ369553.1"="12",
  "OZ369554.1"="X","OZ369555.1"="Y",
  # Anax junius
  "OZ349800.1"="1","OZ349801.1"="2","OZ349802.1"="3","OZ349803.1"="4",
  "OZ349804.1"="5","OZ349805.1"="6","OZ349806.1"="7","OZ349807.1"="8",
  "OZ349808.1"="9","OZ349809.1"="10","OZ349810.1"="11","OZ349811.1"="12",
  "OZ349812.1"="13","OZ349813.1"="X",
  # Archilestes grandis
  "OZ349685.1"="1","OZ349686.1"="2","OZ349687.1"="3","OZ349688.1"="4",
  "OZ349689.1"="5","OZ349690.1"="6","OZ349691.1"="7","OZ349692.1"="8",
  "OZ349693.1"="9","OZ349694.1"="10","OZ349695.1"="11","OZ349696.1"="X",
  # Cordulegaster boltonii
  "OZ297306.1"="1","OZ297307.1"="2","OZ297308.1"="3","OZ297309.1"="4",
  "OZ297310.1"="5","OZ297311.1"="6","OZ297312.1"="7","OZ297313.1"="8",
  "OZ297314.1"="9","OZ297315.1"="10","OZ297316.1"="11","OZ297317.1"="12",
  "OZ297318.1"="X",
  # Ischnura elegans
  "NC_060246.1"="1","NC_060247.1"="2","NC_060248.1"="3","NC_060249.1"="4",
  "NC_060250.1"="5","NC_060251.1"="6","NC_060252.1"="7","NC_060253.1"="8",
  "NC_060254.1"="9","NC_060255.1"="10","NC_060256.1"="11","NC_060257.1"="12",
  "NC_060258.1"="13","NC_060259.1"="X",
  "NW_025791657.1"="?","NW_025791658.1"="?","NW_025791659.1"="?",
  "NW_028141064.1"="?",
  # Orthetrum coerulescens
  "OZ294519.1"="1","OZ294520.1"="2","OZ294521.1"="3","OZ294522.1"="4",
  "OZ294523.1"="5","OZ294524.1"="6","OZ294525.1"="7","OZ294526.1"="8",
  "OZ294527.1"="9","OZ294528.1"="10","OZ294529.1"="11","OZ294530.1"="12",
  "OZ294531.1"="X","OZ294532.1"="Y",
  # Platycnemis pennipes
  "OW121847.1"="1","OW121848.1"="2","OW121849.1"="3","OW121850.1"="4",
  "OW121851.1"="5","OW121852.1"="6","OW121853.1"="7","OW121854.1"="8",
  "OW121855.1"="9","OW121856.1"="10","OW121857.1"="11","OW121858.1"="12",
  "OW121846.1"="X",
  # Pantala flavescens
  "CM036878.1"="1","CM036879.1"="2","CM036880.1"="3","CM036881.1"="4",
  "CM036882.1"="5","CM036883.1"="6","CM036884.1"="7","CM036885.1"="8",
  "CM036886.1"="9","CM036887.1"="10","CM036888.1"="11","CM036889.1"="X",
  # Somatochlora arctica
  "OZ407513.1"="1","OZ407514.1"="2","OZ407515.1"="3","OZ407516.1"="4",
  "OZ407517.1"="5","OZ407518.1"="6","OZ407519.1"="7","OZ407520.1"="8",
  "OZ407521.1"="9","OZ407522.1"="10","OZ407523.1"="11","OZ407524.1"="12",
  "OZ407525.1"="X",
  # Sympetrum striolatum
  "OX388341.1"="1","OX388342.1"="2","OX388343.1"="3","OX388344.1"="4",
  "OX388345.1"="5","OX388346.1"="6","OX388347.1"="7","OX388348.1"="8",
  "OX388349.1"="9","OX388350.1"="10","OX388351.1"="11","OX388352.1"="X"
)

chrLabFun <- function(x) ifelse(x %in% names(map), unname(map[x]), x)
drop_chr  <- names(map)[map %in% c("?", "Y")]   # hide unplaced + Y scaffolds

pal    <- colorRampPalette(c("darkorange","skyblue","darkblue",
                             "purple","darkred","salmon"))
addThm <- ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white"))

# phased blocks for this reference: run plot_riparian once if the csv is absent
blkFile <- file.path(gs$paths$riparian, sprintf("%s_phasedBlks.csv", refGenome))
if (!file.exists(blkFile))
  GENESPACE::plot_riparian(
    gsParam = gs, refGenome = refGenome, genomeIDs = genomeIDs,
    pdfFile = tempfile(fileext = ".pdf")
  )

# drop the hidden scaffolds from the genes and the blocks
bed  <- GENESPACE::read_combBed(file.path(gs$paths$results, "combBed.txt"))
blks <- data.table::fread(blkFile)
bed  <- bed[!(chr  %in% drop_chr)]
blks <- blks[!(chr1 %in% drop_chr | chr2 %in% drop_chr)]

# draw the riparian plot
svglite::svglite("odonate.synteny.svg", width = 11, height = 6.2)
pltDat <- GENESPACE::riparian_engine(
  blk               = blks,
  bed               = bed,
  refGenome         = refGenome,
  genomeIDs         = genomeIDs,
  labelTheseGenomes = genomeIDs,
  palette           = pal,
  chrLabFun         = chrLabFun,
  chrBorderLwd      = 0.1,
  chrExpand         = 2,
  chrFill           = "white",
  howSquare         = 5,
  chrLabFontSize    = 10,
  chrBorderCol      = "black",
  syntenyWeight     = 1,
  addThemes         = addThm
)
print(pltDat$ggplotObj)
dev.off()
