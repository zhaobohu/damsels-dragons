# GENESPACE run for the 29 palaeopteran genomes (inputs in rawGenomeRepo/, built by genespace.sbatch)

library(GENESPACE)
library(data.table)

# make MCScanX failures non-fatal: a divergent pair with zero collinear blocks
# makes MCScanX exit 1; treat that case as "no blocks"
local({
  orig <- getFromNamespace("run_mcscanx", "GENESPACE")
  safe <- function(...) {
    out <- tryCatch(orig(...), error = function(e) {
      if (grepl("External command failed with exit code 1",
                conditionMessage(e), fixed = TRUE))
        return(stats::setNames(integer(0), character(0)))
      stop(e)
    })
    if (is.null(out)) stats::setNames(integer(0), character(0)) else out
  }
  assignInNamespace("run_mcscanx", safe, ns = "GENESPACE")
})

spp <- c(
  "Aeshna_juncea",
  "Somatochlora_meridionalis",
  "Somatochlora_arctica",
  "Serratella_ignita",
  "Anax_junius",
  "Sympetrum_sanguineum",
  "Tanypteryx_hageni",
  "Ischnura_elegans",
  "Cordulegaster_boltonii",
  "Pyrrhosoma_nymphula",
  "Orthetrum_coerulescens",
  "Anax_strenuus",
  "Sympetrum_striolatum",
  "Parafronurus_youi",
  "Ischnura_senegalensis",
  "Platycnemis_pennipes",
  "Ephemera_vulgata",
  "Ephemera_orientalis",
  "Ceriagrion_tenellum",
  "Ecdyonurus_torrentis",
  "Siphlonurus_alternatus",
  "Cloeon_dipterum",
  "Heteragrion_aequatoriale",
  "Baetis_rhodani",
  "Baetis_atlanticus",
  "Archilestes_grandis",
  "Procloeon_bifidum",
  "Hetaerina_titia",
  "Pantala_flavescens"
)

nCores <- as.integer(16)

# gff + proteins -> GENESPACE bed/peptide files
parse_annotations(
  rawGenomeRepo    = "rawGenomeRepo",
  genomeIDs        = spp,
  genomeDirs       = spp,
  genespaceWd      = "genespace",
  presets          = "none",
  gffIdColumn      = "ID",
  headerEntryIndex = 1,
  headerSep        = " ",
  dropDuplicates   = TRUE
)

# set up and run GENESPACE; the state is saved for the riparian figure
gs <- init_genespace(
  wd               = "genespace",
  genomeIDs        = spp,
  path2mcscanx     = Sys.getenv("MCSCANX", "MCScanX"),
  diamondUltraSens = TRUE,
  blkSize          = 3,
  nGaps            = 30,
  onlyOgAnchors    = FALSE,
  nCores           = nCores,
  verbose          = TRUE
)

gs <- run_genespace(gs)
saveRDS(gs, file.path("genespace", "state.rds"))
