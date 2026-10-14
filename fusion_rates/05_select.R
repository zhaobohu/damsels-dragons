#!/usr/bin/env Rscript
# BLAST coverage filter per marker -> <marker>.fna (COI also deduplicated to one record per species)
# pident is not used: identity to the A. juncea reference would bias the tree toward Aeshna

suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
source("config.R")

HITCOLS <- c("qseqid", "qlen", "slen", "qstart", "qend", "sstart", "send", "pident")

# total length covered by a set of intervals, overlaps merged (a record can have several HSPs)
union_len <- function(st, en) {
  a <- pmin(st, en); b <- pmax(st, en)
  o <- order(a); a <- a[o]; b <- b[o]
  tot <- 0L; cs <- a[1]; ce <- b[1]
  if (length(a) > 1) for (i in 2:length(a)) {
    if (a[i] <= ce + 1L) { ce <- max(ce, b[i]) }
    else { tot <- tot + (ce - cs + 1L); cs <- a[i]; ce <- b[i] }
  }
  as.integer(tot + (ce - cs + 1L))
}

# per-record coverage of the reference (scov) and of the record itself (qcov) from hits_<marker>.tsv
read_hits <- function(m) {
  f <- sprintf("hits_%s.tsv", m)
  if (!file.exists(f) || file.size(f) == 0) return(NULL)
  h <- fread(f, col.names = HITCOLS)
  h[, .(scov = union_len(sstart, send) / slen[1],    # of the REFERENCE
        qcov = union_len(qstart, qend) / qlen[1],    # of the RECORD
        alen = union_len(qstart, qend)), by = qseqid]
}

# pass: covers most of the reference, or lies almost wholly inside it (partial barcodes)
pass_cov <- function(d) d$scov > SC_MIN | (d$qcov > QCOV_MIN & d$alen >= MIN_ALN)

kept <- list(); rep <- list()

# COI: filter the pooled records, then keep the best record per accepted species name
if (file.exists("screened.COI.fna")) {
  s <- readDNAStringSet("screened.COI.fna"); names(s) <- sub("^_R_", "", names(s))
  p <- fread("prep.csv")[record_id %chin% names(s)]
  sc <- fread("screen_log.tsv")[, .(record_id, unamb_len, unamb_prop)]
  p  <- merge(p, sc, by = "record_id", all.x = TRUE)
  h  <- read_hits("COI")
  stopifnot(!is.null(h))
  d <- merge(p, h, by.x = "record_id", by.y = "qseqid")
  cat(sprintf("COI: %d of %d screened records had no BLAST hit\n", nrow(p) - nrow(d), nrow(p)))
  stopifnot(nrow(d) > 0.9 * nrow(p))          # catches a truncated / timed-out blast
  # best first, then grouped by acc_name (not the source name, so synonyms share one block)
  d <- d[order(-scov, -unamb_len, -unamb_prop)]
  d <- d[order(family_bold, acc_name)]
  d <- d[pass_cov(d)][!duplicated(acc_name)]
  stopifnot(nrow(d) > 0, !anyDuplicated(d$acc_name))
  sel <- s[d$record_id]; names(sel) <- gsub(" ", "_", d$acc_name)
  writeXStringSet(sel, "COI.fna")
  fwrite(d, "input.spp.csv")
  kept[["COI"]] <- d$acc_name
  rep[[length(rep) + 1]] <- data.table(marker = "COI", queried = uniqueN(p$acc_name),
    with_hit = uniqueN(d$acc_name), kept = nrow(d),
    med_scov = round(median(d$scov), 3),
    via_qcov = sum(d$scov <= SC_MIN))
}

# other markers: one record per species, kept or not
for (m in setdiff(MARKERS, "COI")) {
  f <- sprintf("raw_%s.fasta", m)
  if (!file.exists(f) || file.size(f) == 0) next
  s <- readDNAStringSet(f); names(s) <- sub("^_R_", "", names(s))
  h <- read_hits(m)
  if (is.null(h)) { cat("skip", m, "(no hits file)\n"); next }
  h <- h[qseqid %chin% names(s)]
  k <- h[pass_cov(h)]
  if (!nrow(k)) { cat("skip", m, "(nothing passed coverage)\n"); next }
  writeXStringSet(s[k$qseqid], sprintf("%s.fna", m))
  kept[[m]] <- gsub("_", " ", k$qseqid)
  rep[[length(rep) + 1]] <- data.table(marker = m, queried = length(s),
    with_hit = nrow(h), kept = nrow(k), med_scov = round(median(k$scov), 3),
    via_qcov = sum(k$scov <= SC_MIN))
}

# coverage report
rp <- rbindlist(rep)
rp[, lost_to_coverage := with_hit - kept]
print(rp)
fwrite(rp, "coverage_report.tsv", sep = "\t")

# scov histogram per marker
pdf("scov.pdf", width = 7, height = 5)
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
for (m in MARKERS) {
  h <- read_hits(m); if (is.null(h)) next
  hist(h$scov, breaks = 40, main = m, xlab = "scov", col = "grey80", border = NA)
  abline(v = SC_MIN, col = "red")
}
invisible(dev.off())

# the taxon set the tree will have
tx <- fread(TAXA_TSV)
allsp <- sort(unique(unlist(kept)))
writeLines(gsub(" ", "_", allsp), "tree_taxa.txt")

cat(sprintf("\n%d species have at least one marker after filtering\n", length(allsp)))
cat(sprintf("  with COI          : %d\n", length(kept[["COI"]])))
cat(sprintf("  WITHOUT COI       : %d  (present only via nuclear/other mt markers)\n",
            length(setdiff(allsp, kept[["COI"]]))))
cat(sprintf("karyotyped species represented: %d of %d\n",
            tx[species %chin% allsp, .N], nrow(tx)))
cat(sprintf("  neo-sex among them: %d (%s: %d)\n",
            tx[species %chin% allsp & sex_system != "X", .N], FOCAL_GENUS,
            tx[species %chin% allsp & sex_system != "X" & genus == FOCAL_GENUS, .N]))

# markers per species
n_marks <- list()
for (m in names(kept)) n_marks[[m]] <- data.table(species = kept[[m]], marker = m)
n_marks <- rbindlist(n_marks)[, .N, by = species]
cat("\nmarkers per species:\n"); print(n_marks[, .N, by = N][order(N)])
