# Paths and parameters for the whole pipeline (sourced by every script)

# inputs (read-only, outside the run directory)
KARYO_CSV    <- "/path/to/project/data/odonata.karyotypes.csv"
BOLD_TSV     <- "/path/to/project/data/dragonfly.COI.tsv"
BACKBONE_TRE <- "/path/to/project/data/trees/odonata.FigTree_DivergenceTimes.tre"

# outputs, written into the run directory
TAXA_TSV     <- "taxa.tsv"

# stochastic mapping (09b-09e): maps per dated tree (one map per array task) and RNG offset
NSIM_PER_TREE <- 1
SEED0         <- 20260824

# dating (09): chronos clock model
CHRONOS_MODEL <- "correlated"

BB_TIME_SCALE <- 100

# root IQ-TREE trees on this suborder (sister to Anisoptera + Anisozygoptera)
OUTGROUP_SUBORDER <- "Zygoptera"     # matched against the suborder column of taxa.tsv

# refuse a focal MRCA clade larger than this multiple of the focal tip count
MAX_FOCAL_BLOAT <- 2.0

# focal clade: the 15 Aeshna, excluding A. isoceles (reassigned to Isoaeshna in reconcile.R)
FOCAL_GENUS <- "Aeshna"
FOCAL_TIPS  <- c("Aeshna caerulea", "Aeshna canadensis", "Aeshna clepsydra",
                 "Aeshna crenata", "Aeshna cyanea", "Aeshna grandis",
                 "Aeshna juncea", "Aeshna mixta", "Aeshna palmata",
                 "Aeshna serrata", "Aeshna subarctica", "Aeshna umbrosa",
                 "Aeshna verticalis", "Aeshna viridis", "Aeshna walkeri")

# name resolution (00): orders a GBIF match may land in (the backbone has mayfly outgroups)
ALLOWED_ORDERS <- c("Odonata", "Ephemeroptera")

# markers used for the tree
MARKERS_TREE  <- c("COI", "EF1a")

# COI screening (03): minimum unambiguous (ACGT) length
MIN_UNAMB <- 300

# COI screening (03): keep a record if its longest stop-free ORF (best of six frames)
# covers at least this fraction of it
ORF_FRAC <- 0.80

# references (04): A. juncea mitogenome for the COI Folmer cut
MITO_ASM  <- "/path/to/project/data/genomes/Aeshna.juncea.fna"
MITO_SEQ  <- "OZ369556.1"
PRIMER_MM <- 8            # primer mismatch budget for the Folmer cut

# references (04): non-COI markers use the record at this length quantile (keep it short,
# so short records can still pass the coverage filter)
REF_QUANTILE_DEFAULT <- 0.05
REF_QUANTILE <- list()    # e.g. list(`28S` = 0.25) to override one marker

# coverage filter (05): pass if the record covers SC_MIN of the reference,
# or QCOV_MIN of the record aligned with at least MIN_ALN bases
SC_MIN   <- 0.7
QCOV_MIN <- 0.90
MIN_ALN  <- 280

# sequence retrieval (02): markers downloaded, BLASTed and filtered
MARKERS      <- MARKERS_TREE
# NCBI identification; an API key (export ENTREZ_KEY=<key>) raises the rate limit
ENTREZ_EMAIL <- Sys.getenv("ENTREZ_EMAIL")   # export ENTREZ_EMAIL=<your address> before running
ENTREZ_TOOL  <- "odonata-fusion-pipeline"
ENTREZ_KEY   <- Sys.getenv("ENTREZ_KEY")
# allowed sequence length window per marker (bp)
MARKER_LEN <- list(
  COI  = c(400, 1800),   # barcode fragment up to full-length COI
  `16S`= c(350,  1700),
  `12S`= c(300,  1100),
  `28S`= c(400,  4000),  # D1-D2 through near-full 28S
  `18S`= c(400,  2500),
  EF1a = c(300,  1500),
  H3   = c(250,   500)   # H3 is ~410 bp; anything longer is not the gene
)
MIN_LEN      <- 300                          # fallback if a marker is missing above
N_COI        <- 3           # NCBI COI records per gap species
N_MARKER     <- 3           # records per species for the other markers;
                            # 05_select.R picks among them on coverage
# esearch chunking: term-length budget (NCBI rejects long GET URLs) and species ceiling
MAX_TERM_CHARS <- 2500
BATCH_SIZE     <- 120     # ceiling; the length budget normally bites first

# species per query in the sweep pass that re-queries batch misses (must be 1; 0 skips it)
SWEEP_SIZE     <- 1

# retries for transient NCBI errors
RETRY_TRIES    <- 4      # total attempts per call
RETRY_WAIT     <- 2      # seconds before the first retry, doubling thereafter
MAX_LEN      <- 5000

# tree inference (08, 08b): UFBoot replicates and threads
N_BOOT   <- 2000
THREADS  <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))

# test parameters (13_figures.R): candidate clade size window, in karyotyped tips
N_MIN <- 8
N_MAX <- 24

K_MIN <- 3      # minimum fusion events for a clade to enter the max
SEED  <- 20260824
MC_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

# null probability that a fusion joins the sex chromosome and an autosome, given haploid nchr
# (Schmalz & Blackmon 2026 weighted pSAF; weighted = FALSE gives Anderson et al. 2020)
p_SA <- function(nchr, sys = "XO", weighted = TRUE) {
  l    <- ifelse(sys == "XY" | sys == "neo", 1, 0)    # Y in males
  Da   <- 2 * (nchr - 1)                              # diploid autosomes
  Dhom <- Da + 2                                      # XX females
  Dhet <- Da + 1 + l                                  # X0 or XY males
  SA_hom <- 1 - Da * (Da - 2) / (Dhom * (Dhom - 2))   # one X pair: no X-X fusion
  SA_het <- 1 - Da * (Da - 2) / (Dhet * (Dhet - 2))   # one X, at most one Y
  P_SA <- 0.5 * SA_hom + 0.5 * SA_het
  if (!weighted) return(P_SA)
  P_XA <- 0.5 * SA_hom + 0.5 * SA_het / (1 + l)       # male SA split X:Y as in pSAF
  P_YA <- 0.5 * SA_het * l / (1 + l)
  num  <- 4 / 3 * P_XA + 4 * P_YA
  num / ((1 - P_SA) + num)                            # AA weight 1
}

# null probability used by the tests: SA fusions are scored only in X0 lineages, so neo gets 0
p_null <- function(nchr, sys) ifelse(sys == "XO", p_SA(nchr, sys), 0)
