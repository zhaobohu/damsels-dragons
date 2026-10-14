#!/usr/bin/env Rscript
# Screen the pooled COI records before alignment: unambiguous length and coding
# (longest stop-free ORF over six frames) -> screened.COI.fna, screen_log.tsv

suppressPackageStartupMessages({ library(data.table); library(Biostrings) })
source("config.R")

# records and metadata; duplicated ids mean an upstream step appended twice
s <- readDNAStringSet("unfiltered.COI.fna")
names(s) <- sub("^_R_", "", names(s))
p <- fread("prep.csv")
stopifnot(all(names(s) %in% p$record_id))
if (anyDuplicated(names(s)))
  stop(sprintf("unfiltered.COI.fna has %d duplicated record ids (e.g. %s) -- rerun 01_prep.R then 02_fetch.R",
               sum(duplicated(names(s))), names(s)[duplicated(names(s))][1]))

# unambiguous (ACGT) length
cln <- DNAStringSet(gsub("[^ACGTNacgtn]", "", as.character(s)))
names(cln) <- names(s)
unamb <- rowSums(letterFrequency(cln, c("A", "C", "G", "T")))
prop  <- unamb / width(cln)

gc5 <- getGeneticCode("5")     # invertebrate mitochondrial

# longest stop-free codon run (nt) in frame f
orf_nt <- function(x, f) {
  w <- (width(x) - f) %/% 3L * 3L
  keep <- w > 0
  out <- rep(0L, length(x))
  if (!any(keep)) return(out)
  aa <- as.character(translate(subseq(x[keep], start = f + 1L, width = w[keep]),
                               genetic.code = gc5, if.fuzzy.codon = "solve"))
  out[keep] <- vapply(strsplit(aa, "*", fixed = TRUE),
                      function(p) if (!length(p)) 0L else max(nchar(p)), 0L) * 3L
  out
}

# internal stop codons in frame f (for the report only)
n_stop <- function(x, f) {
  w <- (width(x) - f) %/% 3L * 3L
  keep <- w > 0
  out <- rep(NA_integer_, length(x))
  if (!any(keep)) return(out)
  aa <- sub("[*]+$", "",
        as.character(translate(subseq(x[keep], start = f + 1L, width = w[keep]),
                               genetic.code = gc5, if.fuzzy.codon = "solve")))
  out[keep] <- vapply(gregexpr("*", aa, fixed = TRUE),
                      function(g) if (g[1] == -1L) 0L else length(g), 0L)
  out
}

# fewest stops over the six frames (3 forward, 3 reverse complement)
rc <- reverseComplement(cln)
stops6 <- cbind(vapply(0:2, function(f) n_stop(cln, f), integer(length(cln))),
                vapply(0:2, function(f) n_stop(rc,  f), integer(length(cln))))
min_stops <- apply(stops6, 1, function(v)
  if (all(is.na(v))) NA_integer_ else min(v, na.rm = TRUE))

# the filter: longest stop-free ORF over the six frames as a fraction of the record;
# records that pass only on the reverse strand are reverse-complemented
orfs6 <- cbind(vapply(0:2, function(f) orf_nt(cln, f), integer(length(cln))),
               vapply(0:2, function(f) orf_nt(rc,  f), integer(length(cln))))
best_orf <- apply(orfs6, 1, max)
best_i   <- apply(orfs6, 1, which.max)
orf_frac <- best_orf / pmax(width(cln), 1L)

coding_ok <- orf_frac >= ORF_FRAC
flip <- coding_ok & best_i > 3L
if (any(flip)) cln[flip] <- rc[flip]
stopfree <- coding_ok            # stop_free column in the log

# log every record with the reason it was dropped, write the kept ones
keep <- unamb > MIN_UNAMB & stopfree
log <- data.table(record_id = names(s), unamb_len = unamb, unamb_prop = round(prop, 3),
                  min_stops = min_stops, stop_free = stopfree,
                  strand = fifelse(flip, "-", "+"), kept = keep)
log[, reason := fifelse(keep, "", fifelse(!stopfree & unamb <= MIN_UNAMB,
      "short+numt", fifelse(!stopfree, "no stop-free frame", "too short")))]
fwrite(log, "screen_log.tsv", sep = "\t")

writeXStringSet(cln[keep], "screened.COI.fna")

# report
cat(sprintf("%d records in\n", length(s)))
cat(sprintf("  dropped %d with unambiguous length <= %d\n", sum(unamb <= MIN_UNAMB), MIN_UNAMB))
cat(sprintf("  dropped %d whose longest stop-free ORF covers < %.0f%% of the record\n",
            sum(!coding_ok), 100 * ORF_FRAC))
cat(sprintf("  reverse-complemented %d minus-strand records (kept, not dropped)\n", sum(flip)))

# ORF fractions and stop counts of the coding rejects
sel <- !coding_ok & unamb > MIN_UNAMB
if (any(sel)) {
  cat(sprintf("\nORF fraction among the %d rejected:\n", sum(sel)))
  print(table(cut(orf_frac[sel], c(0, .3, .5, .6, .7, .8), include.lowest = TRUE), dnn = NULL))
  cat(sprintf("  their stop counts: median %d, range %d-%d\n",
              median(min_stops[sel]), min(min_stops[sel]), max(min_stops[sel])))
  cat("  Lower ORF_FRAC in config.R only if the top bin looks like real COI.\n")
}
cat(sprintf("%d records kept, covering %d species\n", sum(keep),
            p[record_id %in% names(s)[keep], uniqueN(acc_name)]))

# species that lost every record, and why
lost <- setdiff(p[record_id %in% names(s), unique(acc_name)],
                p[record_id %in% names(s)[keep], unique(acc_name)])
if (length(lost)) {
  cat(sprintf("\n%d species lost every record:\n", length(lost)))
  li <- data.table(record_id = names(s), acc = p$acc_name[match(names(s), p$record_id)],
                   unamb = unamb, ms = min_stops, of = round(orf_frac, 2))[acc %chin% lost]
  print(li[, .(n_records = .N, max_unamb_len = max(unamb),
               min_stops = min(ms, na.rm = TRUE), best_orf_frac = max(of)), by = acc][order(acc)],
        nrows = 30)
  cat("max_unamb_len <= MIN_UNAMB means it was dropped for length, not coding.\n")
}
