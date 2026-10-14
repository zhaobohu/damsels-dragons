library(data.table)

run_dir <- "../gene_pairs_ds"

source("../lib/anchor_filter.R")
min.ratio     <- 0.65   # keep equal to figs/fig3/fig3c.R
min.pid       <- 0
ratio.na.keep <- FALSE
infile <- file.path(run_dir, "all_ds_results.tsv")
g_min <- 0.5; f_min <- 0.3
br <- as.numeric(readLines(file.path(run_dir, "breakpoint.txt")))

pg <- fread("../odonate_synteny/genespace/pangenes/Aeshna_juncea_pangenes.txt.gz")

# strip isoform suffix so IDs match across GFF / pangenome / dS table
# (A. juncea IDs carry a legacy "anax_" prefix; species always comes from the pangenome's genome column)
gid <- function(x) sub("\\.t[0-9]+$", "", x)

# dS table, with the QC columns as numbers
DT <- fread(infile)
nc <- c("dS", "mean_guidance_score", "frac_cols_retained")
for (col in nc) {
  DT[, (col) := suppressWarnings(as.numeric(get(col)))]
}

# reciprocal 1:1 MCScanX anchors of one comparison that pass the rank filter (anchor_filter.R)
recip_anchors <- function(cmp) {
  parts <- strsplit(cmp, "__vs__")[[1]]
  p   <- file.path(run_dir, "comp", cmp, "mcscanx.pairs.tsv")
  ida <- file.path(run_dir, "chrom", paste0(parts[1], ".ids"))
  if (!file.exists(p) || file.size(p) == 0 || !file.exists(ida)) {
    warning("no MCScanX anchors for ", cmp, " (", p, ")")
    return(NULL)
  }
  a   <- fread(p, header = FALSE, select = 1:2); setnames(a, c("g1", "g2"))
  ids <- readLines(ida)
  s1  <- a$g1 %in% ids
  # transcript-level IDs until after the rank filter (ab.rank.blast is keyed on them)
  d   <- unique(data.table(gene_a = fifelse(s1, a$g1, a$g2),
                           gene_b = fifelse(s1, a$g2, a$g1)))
  d[, `:=`(ka = gid(gene_a), kb = gid(gene_b))]
  d   <- d[, na := .N, by = ka][, nb := .N, by = kb][na == 1L & nb == 1L]
  if (!nrow(d)) return(NULL)

  n0 <- nrow(d)
  d  <- attach_rank(d, cmp, run_dir)
  cat(sprintf("[anchors] %s: %d -> %d after rank filter (min.ratio = %.2f)\n",
              cmp, n0, nrow(d), min.ratio))
  rank_report(d, cmp)
  unique(data.table(gene_a = gid(d$gene_a), gene_b = gid(d$gene_b)))
}

# Step 1: QC-filtered X-Y pairs that are reciprocal 1:1 anchors
xy_qc <- DT[comparison == "juncea_X__vs__juncea_Y"][
  is.finite(dS) & is.finite(mean_guidance_score) & mean_guidance_score >= g_min &
  is.finite(frac_cols_retained) & frac_cols_retained >= f_min]
cat(sprintf("X-Y pairs after QC (guidance>=%.2f, cols>=%.2f): %d\n",
            g_min, f_min, nrow(xy_qc)))

anchors <- recip_anchors("juncea_X__vs__juncea_Y")
xy <- xy_qc[, .(g_x = gid(gene_a), g_y = gid(gene_b), dS)]
gametologs <- xy[anchors, on = .(g_x = gene_a, g_y = gene_b), nomatch = 0]
cat(sprintf("X-Y pairs after reciprocal 1:1 restriction: %d\n\n", nrow(gametologs)))

# Step 2: keep gametologs whose X gene is on the neo-X arm (OZ369554.1, start < breakpoint)
pg_ajun <- pg[genome == "Aeshna_juncea"]
pg_ajun[, id_stripped := gid(id)]

neox_x_genes <- pg_ajun[chr == "OZ369554.1" & start < br, id_stripped]
gametologs_neox <- gametologs[g_x %in% neox_x_genes]
cat(sprintf("  %d pairs have X gene on neo-X arm (OZ369554.1, start < %s bp)\n",
            nrow(gametologs_neox), format(br, big.mark = ",")))

# Step 3: pangene ID of each X gene
gametologs_neox[, pg_x := pg_ajun[match(g_x, id_stripped), pgID]]
n_missing <- sum(is.na(gametologs_neox$pg_x))
if (n_missing > 0) {
  cat(sprintf("  %d X genes not in pangenome (dropped)\n", n_missing))
  gametologs_neox <- gametologs_neox[!is.na(pg_x)]
}
cat(sprintf("  %d candidate pangenes on neo-X arm\n", nrow(gametologs_neox)))

# Step 4: outgroup orthologs per candidate pangene
outgroups <- c("Anax_junius", "Somatochlora_arctica", "Cordulegaster_boltonii",
               "Sympetrum_striolatum", "Orthetrum_coerulescens", "Pantala_flavescens",
               "Ischnura_elegans", "Platycnemis_pennipes", "Archilestes_grandis")

candidate_pg_data <- pg[pgID %in% gametologs_neox$pg_x]
if (nrow(candidate_pg_data) == 0) stop("no candidate pangenes — check upstream ID matching")

counts <- dcast(candidate_pg_data, pgID ~ genome,
                fun.aggregate = length, value.var = "id")

# every outgroup gets a column, even with no genes in any candidate pangene
for (sp in outgroups) if (!sp %in% names(counts)) counts[[sp]] <- 0L

# Step 5: keep pangenes with exactly 1 ortholog in every outgroup
one_each <- rep(TRUE, nrow(counts))
for (sp in outgroups) {
  one_each <- one_each & counts[[sp]] == 1
}
good_pgs <- counts[one_each, pgID]
locus_set <- gametologs_neox[pg_x %in% good_pgs]
cat(sprintf("  %d pangenes have 1:1 orthologs across all %d outgroups\n",
            nrow(locus_set), length(outgroups)))

# per-outgroup retention
cat("\nPer-outgroup retention (# pangenes with exactly 1 ortholog):\n")
for (sp in outgroups) {
  cat(sprintf("  %-30s %d / %d\n", sp,
              sum(counts[[sp]] == 1), nrow(counts)))
}

# Step 6: 11-taxon gene ID table, one row per locus
outgroup_genes <- dcast(
  candidate_pg_data[pgID %in% good_pgs & genome %in% outgroups],
  pgID ~ genome, value.var = "id"
)

locus_table <- merge(locus_set[, .(pg_x, g_x, g_y, dS)],
                     outgroup_genes,
                     by.x = "pg_x", by.y = "pgID")

cat(sprintf("\nFinal locus table: %d loci x 11 tips\n", nrow(locus_table)))
fwrite(locus_table, "locus_table.tsv", sep = "\t")
