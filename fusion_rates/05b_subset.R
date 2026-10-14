#!/usr/bin/env Rscript
# Restrict to MARKERS_TREE before 06: stash other markers in subset_stash/, rewrite tree_taxa.txt

suppressPackageStartupMessages({ library(Biostrings); library(data.table) })
source("config.R")

stopifnot(file.exists(TAXA_TSV))
tx <- fread(TAXA_TSV)
tx[, tip := gsub(" ", "_", species)]

# per-marker files written by 05_select.R
fna <- setdiff(list.files(".", "\\.fna$"), c("unfiltered.COI.fna", "screened.COI.fna"))
if (!length(fna)) stop("no <marker>.fna here -- run 05_select.R first", call. = FALSE)

keep_m <- MARKERS_TREE
drop_m <- setdiff(sub("\\.fna$", "", fna), keep_m)

# taxa allowed in the subset (all of taxa.tsv)
tips_ok <- tx$tip

# move dropped markers (and their trimmed alignments) to subset_stash/
dir.create("subset_stash", showWarnings = FALSE)
for (m in drop_m) {
  file.rename(sprintf("%s.fna", m), sprintf("subset_stash/%s.fna", m))
  a <- sprintf("align/%s.trim.fasta", m)
  if (file.exists(a)) file.rename(a, sprintf("subset_stash/%s.trim.fasta", m))
}
if (length(drop_m)) cat(sprintf("stashed %d marker(s): %s\n",
                                length(drop_m), paste(drop_m, collapse = ", ")))

# kept markers: restrict to allowed tips and delete their now-stale trimmed alignments
present <- character(0)
for (m in keep_m) {
  f <- sprintf("%s.fna", m)
  if (!file.exists(f)) { cat(sprintf("  %s.fna not found, skipping\n", m)); next }
  s <- readDNAStringSet(f)
  names(s) <- sub("^_R_", "", names(s))
  k <- s[names(s) %chin% tips_ok]
  if (!length(k)) { cat(sprintf("  %s: no tips survive the taxon filter\n", m)); next }
  writeXStringSet(k, f)
  a <- sprintf("align/%s.trim.fasta", m)
  if (file.exists(a)) file.remove(a)
  present <- union(present, names(k))
  cat(sprintf("  %-5s %3d -> %3d tips\n", m, length(s), length(k)))
}
if (!length(present)) stop("nothing left after subsetting -- check MARKERS_TREE",
                           call. = FALSE)

writeLines(sort(present), "tree_taxa.txt")

# report the retained taxon set
nb <- tx[tip %chin% present]
cat(sprintf("\n%d tips retained\n", length(present)))
cat(sprintf("  neo-sex          : %d\n", nb[sex_system != "X", .N]))
cat(sprintf("  %s%s: %d (%d neo-sex)\n", FOCAL_GENUS, strrep(" ", max(0, 17 - nchar(FOCAL_GENUS))),
            nb[genus == FOCAL_GENUS, .N], nb[genus == FOCAL_GENUS & sex_system != "X", .N]))
cat(sprintf("  n range          : %d-%d\n", min(nb$nchr, na.rm = TRUE), max(nb$nchr, na.rm = TRUE)))
if (nb[genus == FOCAL_GENUS & sex_system != "X", .N] < 5)
  cat("\nWARNING: fewer than 5 neo-sex ", FOCAL_GENUS, " tips -- that IS the signal.\n", sep = "")
cat("\nalign/ entries for the retained markers were deleted -- rerun 06 onward:\n")
cat("  bash 06_align_trim.sh && Rscript 07_constraint.R && Rscript 07b_check.R\n")

# list outputs built from the previous taxon set that a rerun would not overwrite
# (09_date.R resumes from tree/boot_dated_ckpt.rds; the best-score grep reads every run.*/)
stale <- c(
  "tree/boot_dated_ckpt.rds", "tree/odon.dated.tre", "tree/odon.boot.dated.tre",
  "BEST_TREE", "prep.rds", "simmaps.rds", "simmap_ckpt.rds", "todo.txt",
  "tree/constraint.tre", "tree/partitions.nex")
stale <- stale[file.exists(stale)]
dirs  <- c(list.files(".", "^run\\.[0-9]+$", include.dirs = TRUE),
           if (dir.exists("support")) "support", if (dir.exists("maps")) "maps")
dirs  <- dirs[dir.exists(dirs)]

if (length(stale) || length(dirs)) {
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("STALE ARTIFACTS from the previous taxon set -- DELETE BEFORE RERUNNING\n")
  cat(strrep("=", 70), "\n", sep = "")
  if (length(stale)) cat("  files: ", paste(stale, collapse = ", "), "\n", sep = "")
  if (length(dirs))  cat("  dirs : ", paste(dirs, collapse = ", "), "\n", sep = "")
  cat("\n  rm -rf ", paste(c(stale, dirs), collapse = " "), "\n", sep = "")
  cat("\n  tree/boot_dated_ckpt.rds is the dangerous one: 09_date.R RESUMES\n")
  cat("  from it, silently mixing dated trees from the old taxon set.\n")
  cat("  Stale run.*/ dirs are next: the best-score grep reads all of them.\n")
}
