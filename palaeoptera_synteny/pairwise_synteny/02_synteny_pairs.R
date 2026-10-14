# Pairwise chromosome-level synteny metrics (ARI, nesting, split90) for every pair of genomes.
# Rscript 02_synteny_pairs.R o2o        -> pairs_o2o.tsv (primary)
# Rscript 02_synteny_pairs.R universal  -> pairs_universal.tsv (robustness)

library(data.table)
library(ape)

MARKER_SET <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(MARKER_SET) || !MARKER_SET %in% c("o2o", "universal"))
  stop("usage: Rscript 02_synteny_pairs.R o2o | universal")

# paths - edit here; outputs land in WORKDIR
WORKDIR <- "/path/to/project/palaeoptera_synteny/pairwise_synteny"
GENESPACE <- "/path/to/project/palaeoptera_synteny/genespace_run/genespace"
OFD   <- file.path(GENESPACE, "orthofinder/Results_Aug18")
BED   <- file.path(GENESPACE, "bed")
CHROM <- "chromosomes.tsv"                     # from 00_chromosomes.py
TREE  <- "tree/iqtree/tree.treefile"           # from 01_build_tree.py
OUT   <- paste0("pairs_", MARKER_SET, ".tsv")
minMarkers <- 50

setwd(WORKDIR)
ogd  <- file.path(OFD, "Orthogroups")
orth <- file.path(OFD, "Orthologues")

eph <- c("Baetis_atlanticus", "Baetis_rhodani", "Cloeon_dipterum", "Procloeon_bifidum",
         "Ecdyonurus_torrentis", "Parafronurus_youi", "Ephemera_orientalis",
         "Ephemera_vulgata", "Serratella_ignita", "Siphlonurus_alternatus")
ani <- c("Aeshna_juncea", "Anax_junius", "Anax_strenuus", "Cordulegaster_boltonii",
         "Orthetrum_coerulescens", "Pantala_flavescens", "Somatochlora_arctica",
         "Somatochlora_meridionalis", "Sympetrum_sanguineum", "Sympetrum_striolatum",
         "Tanypteryx_hageni")

for (p in c(CHROM, TREE)) if (!file.exists(p)) stop("missing input: ", p)
for (p in c(BED, if (MARKER_SET == "o2o") orth else ogd))
  if (!dir.exists(p)) stop("missing input directory: ", p)

# 1. gene -> chromosome, restricted to the chromosome set
chrom <- fread(CHROM)[, .(genome, chr)]
bedf  <- list.files(BED, pattern = "\\.bed$", full.names = TRUE)
if (!length(bedf)) stop("no .bed files in ", BED)
bed <- list()
for (f in bedf) {
  bed[[f]] <- data.table(genome = sub("\\.bed$", "", basename(f)),
                         fread(f, header = FALSE, select = c(1, 4), col.names = c("chr", "gene")))
}
bed <- rbindlist(bed)
g2c <- merge(bed, chrom, by = c("genome", "chr"))[, .(genome, gene, chr)]
setkey(g2c, genome, gene)

sp <- sort(unique(bed$genome))
if (!setequal(sp, unique(chrom$genome)))
  stop("chromosomes.tsv and bed/ do not cover the same genomes")
zyg <- setdiff(sp, c(eph, ani))
cat("[genomes] ", length(sp), "   Anisoptera ", length(ani), "  Zygoptera ", length(zyg),
    "  Ephemeroptera ", length(eph), "\n", sep = "")

# 2. universal marker set: single-copy orthogroups, split by genome
if (MARKER_SET == "universal") {
  sco <- trimws(readLines(file.path(ogd, "Orthogroups_SingleCopyOrthologues.txt")))
  sco <- sco[nzchar(sco)]
  og  <- fread(file.path(ogd, "Orthogroups.tsv"), sep = "\t", header = TRUE, na.strings = "")
  setnames(og, 1, "og")
  og  <- og[og %in% sco]
  m   <- melt(og, id.vars = "og", variable.name = "genome", value.name = "gene",
              variable.factor = FALSE)[!is.na(gene)]
  m[, gene := trimws(gene)]
  stopifnot(!any(grepl(",", m$gene)))
  d <- merge(m, g2c, by = c("genome", "gene"))
  L <- split(d[, .(og, chr)], d$genome)
  cat("[markers] universal single-copy: ", length(sco), " orthogroups\n", sep = "")
} else {
  cat("[markers] per-pair one-to-one orthologues\n")
}

# 3. pairwise metrics, for pairs sharing at least minMarkers markers
pr <- as.data.table(t(combn(sp, 2))); setnames(pr, c("a", "b"))
res <- list()
for (i in seq_len(nrow(pr))) {
  a <- pr$a[i]; b <- pr$b[i]

  # shared markers with their chromosome in each genome
  if (MARKER_SET == "universal") {
    x <- merge(L[[a]], L[[b]], by = "og", suffixes = c(".a", ".b"))
    n_raw <- NA_integer_
  } else {
    # one-to-one rows of OrthoFinder's Orthologues_<A>/<A>__v__<B>.tsv (either direction)
    f <- file.path(orth, paste0("Orthologues_", a), paste0(a, "__v__", b, ".tsv"))
    if (!file.exists(f))
      f <- file.path(orth, paste0("Orthologues_", b), paste0(b, "__v__", a, ".tsv"))
    if (!file.exists(f)) stop("no orthologue file for ", a, " x ", b)
    o <- fread(f, sep = "\t", header = TRUE)
    if (!all(c(a, b) %in% names(o)))
      stop("unexpected columns in ", basename(f), ": ", paste(names(o), collapse = ", "))
    o <- o[!grepl(",", get(a)) & !grepl(",", get(b)), .(ga = trimws(get(a)), gb = trimws(get(b)))]
    x <- merge(o, g2c[genome == a, .(ga = gene, chr.a = chr)], by = "ga")
    x <- merge(x, g2c[genome == b, .(gb = gene, chr.b = chr)], by = "gb")
    n_raw <- nrow(o)
  }
  if (nrow(x) < minMarkers) next

  # share of markers on each chromosome's modal partner, in both directions
  tt <- x[, .N, by = .(chr.a, chr.b)]
  ra <- tt[, .(top = max(N), tot = sum(N)), by = chr.a]
  rb <- tt[, .(top = max(N), tot = sum(N)), by = chr.b]
  pa <- ra[, sum(top) / sum(tot)]; pb <- rb[, sum(top) / sum(tot)]
  fav <- if (pa >= pb) ra else rb

  # adjusted Rand index between the two chromosome partitions
  ct <- table(x$chr.a, x$chr.b)
  s  <- sum(choose(ct, 2))
  sa <- sum(choose(rowSums(ct), 2)); sb <- sum(choose(colSums(ct), 2))
  e  <- sa * sb / choose(nrow(x), 2)

  # nest = share on the modal partner in the more favourable direction (displacement = 100 * (1 - nest));
  # split90 = share of chromosomes in that direction with < 90% of their markers on one partner
  res[[length(res) + 1]] <- data.table(a = a, b = b, n = nrow(x), n_raw = n_raw,
                                       na = uniqueN(x$chr.a), nb = uniqueN(x$chr.b),
                                       ari = (s - e) / ((sa + sb) / 2 - e),
                                       nest = max(pa, pb),
                                       split90 = fav[, mean(top / tot < 0.9)])
}
res <- rbindlist(res)
cat("[pairs] ", nrow(res), " of ", nrow(pr), " computed\n", sep = "")
if (MARKER_SET == "o2o") {
  um <- 1 - res$n / res$n_raw
  cat(sprintf("[unmatched] one-to-one genes off the chromosome set: median %.1f%%, max %.1f%%\n",
              100 * median(um), 100 * max(um)))
  if (max(um) > 0.15) warning("> 15% unmatched in at least one pair; check gene IDs")
}

# 4. tree distance and pair class
tr <- read.tree(TREE)
if (is.null(tr$edge.length)) stop("tree has no branch lengths")
miss <- setdiff(sp, tr$tip.label)
if (length(miss)) stop("genomes absent from the tree: ", paste(miss, collapse = ", "))
D <- cophenetic(tr)
res[, dist := D[cbind(a, b)]]

cls  <- function(s) fifelse(s %in% eph, "Eph", "Odo")
subc <- function(s) fifelse(s %in% eph, "Eph", fifelse(s %in% ani, "Aniso", "Zygo"))
res[, class := { u <- cls(a); v <- cls(b)
  fifelse(u == v, paste(u, v, sep = "-"), "Odo-Eph") }]
res[, sub := { u <- subc(a); v <- subc(b)
  fifelse(u == "Eph" & v == "Eph", "Eph-Eph",
  fifelse(u == "Eph" | v == "Eph", "Odo-Eph",
  fifelse(u == v, paste0("within-", u), "Aniso-Zygo"))) }]

# write the table and print a summary by class
setcolorder(res, c("a", "b", "n", "n_raw", "na", "nb", "ari", "nest", "split90",
                   "dist", "class", "sub"))
fwrite(res, OUT, sep = "\t")

cat("\n[by class]\n")
print(res[, .(pairs = .N, markers_med = as.numeric(median(n)),
              dist_min = round(min(dist), 3), dist_max = round(max(dist), 3),
              ari_med = round(median(ari), 3),
              displaced_med = round(100 * (1 - median(nest)), 2)),
          by = class][order(class)])
cat("[wrote] ", OUT, "\n", sep = "")
