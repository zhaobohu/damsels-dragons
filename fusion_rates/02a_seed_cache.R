#!/usr/bin/env Rscript
# One-time: import NCBI sequences fetched before caching existed into ncbi_cache.csv

suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
source("config.R")

rows <- list()

# non-COI markers: raw_<marker>.fasta, headers are species names with underscores
for (m in setdiff(MARKERS, "COI")) {
  f <- sprintf("raw_%s.fasta", m)
  if (!file.exists(f) || file.size(f) == 0) next
  s <- readDNAStringSet(f)
  rows[[length(rows) + 1]] <- data.table(
    species = gsub("_", " ", names(s)), marker = m,
    accession = NA_character_, seq = as.character(s))
  cat(sprintf("%-5s %d sequences\n", m, length(s)))
}

# COI: the records prep.csv marks as source == "ncbi"
if (file.exists("prep.csv") && file.exists("unfiltered.COI.fna")) {
  p <- fread("prep.csv")
  if ("source" %in% names(p) && p[source == "ncbi", .N]) {
    s <- readDNAStringSet("unfiltered.COI.fna")
    n <- p[source == "ncbi"][record_id %chin% names(s)]
    rows[[length(rows) + 1]] <- data.table(
      species = n$acc_name, marker = "COI",
      accession = as.character(n$insdc_acs),
      seq = as.character(s[n$record_id]))
    cat(sprintf("%-5s %d sequences\n", "COI", nrow(n)))
  }
}

if (!length(rows)) stop("nothing to harvest -- run 02_fetch.R at least once first")

# write the cache
ca <- unique(rbindlist(rows))
ca[, len := nchar(seq)]
ca <- ca[len > 0]
fwrite(ca, "ncbi_cache.csv")

# length summary per marker
cat(sprintf("\ncached %d sequences for %d species\n", nrow(ca), uniqueN(ca$species)))
cat("\nlength distribution per marker (check against MARKER_LEN in config.R):\n")
print(ca[, .(n = .N, min = min(len), med = as.integer(median(len)), max = max(len)),
         by = marker][order(marker)])

# records outside the current MARKER_LEN window (02_fetch.R will re-query them)
cat("\nrecords that would now FAIL the MARKER_LEN window and be re-queried:\n")
ca[, ok := {
  r <- MARKER_LEN[[marker[1]]]
  if (is.null(r)) r <- c(MIN_LEN, MAX_LEN)
  len >= r[1] & len <= r[2]
}, by = marker]
print(ca[ok == FALSE, .(n_records = .N, n_species = uniqueN(species)), by = marker])
