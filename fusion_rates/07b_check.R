#!/usr/bin/env Rscript
# Check the IQ-TREE inputs (constraint, alignments, partitions, focal clade, families)
# before submitting 08_tree.sbatch; exits non-zero on any problem

suppressPackageStartupMessages({ library(data.table); library(ape); library(Biostrings) })
source("config.R")
source("tree_utils.R")
source("root_tree.R")

# result printers; bad() counts failures
fail <- 0L
bad  <- function(...) { cat("  FAIL ", sprintf(...), "\n", sep = ""); fail <<- fail + 1L }
ok   <- function(...) cat("  ok   ", sprintf(...), "\n", sep = "")
note <- function(...) cat("  note ", sprintf(...), "\n", sep = "")

# 1. constraint parses, no duplicated tips
cat("constraint\n")
if (!file.exists("tree/constraint.tre")) {
  bad("tree/constraint.tre missing -- run 07_constraint.R")
} else {
  ct <- tryCatch(read.tree("tree/constraint.tre"), error = function(e) e)
  if (inherits(ct, "error")) bad("will not parse: %s", conditionMessage(ct))
  else {
    ok("parses, %d tips, %d internal nodes", Ntip(ct), Nnode(ct))
    if (anyDuplicated(ct$tip.label)) bad("%d duplicated tip labels", sum(duplicated(ct$tip.label)))
  }
}

# 2. alignments: no _R_ prefixes, no duplicated labels, equal widths
cat("\nalignments\n")
af <- list.files("align", "\\.trim\\.fasta$", full.names = TRUE)
if (!length(af)) bad("no align/*.trim.fasta -- run 06_align_trim.sh")
labs <- list()
for (p in af) {
  m <- sub("\\.trim\\.fasta$", "", basename(p))
  s <- readBStringSet(p)
  n <- sub("^_R_", "", names(s))
  if (any(grepl("^_R_", names(s)))) bad("%s: %d labels still carry the MAFFT _R_ prefix",
                                        m, sum(grepl("^_R_", names(s))))
  if (anyDuplicated(n)) bad("%s: %d duplicated labels", m, sum(duplicated(n)))
  if (length(unique(width(s))) != 1L) bad("%s: not aligned -- %d distinct widths",
                                          m, length(unique(width(s))))
  labs[[m]] <- n
}
allsp <- sort(unique(unlist(labs)))
ok("%d markers, %d distinct tips across them", length(labs), length(allsp))

# every alignment tip is in tree_taxa.txt
if (file.exists("tree_taxa.txt")) {
  tt <- sub("^_R_", "", readLines("tree_taxa.txt"))
  ex <- setdiff(allsp, tt)
  if (length(ex)) bad("%d alignment tips absent from tree_taxa.txt (e.g. %s)",
                      length(ex), paste(head(ex, 3), collapse = ", "))
  else ok("every alignment tip is in tree_taxa.txt")
} else bad("tree_taxa.txt missing -- run 05_select.R")

# every constraint tip has sequence; tips outside the constraint float free
if (exists("ct") && !inherits(ct, "error")) {
  miss <- setdiff(ct$tip.label, allsp)
  if (length(miss))
    bad("%d constraint tips have NO sequence (IQ-TREE rejects this): %s",
        length(miss), paste(head(miss, 3), collapse = ", "))
  else ok("every constraint tip has sequence")
  free <- setdiff(allsp, ct$tip.label)
  if (length(free)) note("%d tips are outside the constraint and will float free", length(free))
}

# 3. partition charsets point at existing files and fit inside them
cat("\npartitions\n")
if (!file.exists("tree/partitions.nex")) {
  bad("tree/partitions.nex missing")
} else {
  pl <- grep("charset", readLines("tree/partitions.nex"), value = TRUE)
  ok("%d charsets", length(pl))
  files <- unique(sub(".*= *([^:]+):.*", "\\1", pl))
  for (f in files) {
    if (!file.exists(f)) { bad("charset points at a missing file: %s", f); next }
    w <- width(readBStringSet(f))[1]
    ends <- as.integer(sub(".*[-:] *([0-9]+)\\\\?3?;.*", "\\1",
                           grep(f, pl, fixed = TRUE, value = TRUE)))
    if (any(ends > w, na.rm = TRUE))
      bad("%s is %d cols but a charset ends at %d", basename(f), w, max(ends))
  }
  if (!fail) ok("all charset ranges fit inside their alignments")
}

# 4. focal clade: tips with sequence, neo-sex tips, monophyly in the constraint
cat("\nfocal clade\n")
tx <- fread(TAXA_TSV)
ft <- gsub(" ", "_", FOCAL_TIPS)
have <- intersect(ft, allsp)
cat(sprintf("  %s: %d karyotyped, %d with sequence\n", FOCAL_GENUS, length(ft), length(have)))
neo <- tx[species %chin% gsub("_", " ", allsp) & sex_system != "X"]
cat(sprintf("  neo-sex tips in the matrix: %d (%s: %d)\n", nrow(neo), FOCAL_GENUS,
            neo[genus == FOCAL_GENUS, .N]))
if (neo[genus == FOCAL_GENUS, .N] < 5)
  bad("fewer than 5 neo-sex %s tips -- that IS the signal, stop and check", FOCAL_GENUS)
if (exists("ct") && !inherits(ct, "error")) {
  inct <- intersect(ft, ct$tip.label)
  if (length(inct) >= 2 && !is.monophyletic(ct, inct))
    note("%s is not constrained as a clade (fine -- the data will place it)", FOCAL_GENUS)
}

# focal MRCA clade in every available tree (support tree and each run.N), rooted first
rd <- list.dirs(".", recursive = FALSE, full.names = FALSE)
rd <- rd[grepl("^run\\.[0-9]+$", rd)]
tf <- c("support/support.treefile", sprintf("%s/%s.treefile", rd, rd))
tf <- tf[file.exists(tf)]
if (length(tf)) {
  res <- list()
  for (p in tf) {
    tr <- read.tree(p); tr$tip.label <- sub("^_R_", "", tr$tip.label)
    if (!is.rooted(tr)) tr <- suppressWarnings(root_tree(tr, tx))
    f2 <- intersect(ft, tr$tip.label)
    if (length(f2) < 2 || !is.rooted(tr)) next
    ins <- tr$tip.label[tips_by_node(tr)[[getMRCA(tr, f2)]]]
    res[[length(res) + 1]] <- list(tree = basename(p), n_focal = length(f2),
                                   n_clade = length(ins), extra = setdiff(ins, f2))
  }
  if (length(res)) {
    sz <- vapply(res, function(x) x$n_clade, 0L)
    nf <- res[[1]]$n_focal
    cat(sprintf("  %s across %d tree(s): %d focal tips, MRCA clade %d-%d (bloat %.2f-%.2f)\n",
                FOCAL_GENUS, length(res), nf, min(sz), max(sz), min(sz)/nf, max(sz)/nf))
    # non-focal tips inside the clade, and how much the clade size varies across trees
    ex <- sort(unique(unlist(lapply(res, function(x) x$extra))))
    if (length(ex)) {
      n_in <- vapply(ex, function(e) sum(vapply(res, function(x) e %in% x$extra, TRUE)), 0L)
      note("also inside, in n of %d trees: %s", length(res),
           paste(sprintf("%s (%d)", ex, n_in), collapse = ", "))
      span <- max(sz) / min(sz)
      if (max(sz) > MAX_FOCAL_BLOAT * nf)
        note("membership is UNSTABLE and exceeds MAX_FOCAL_BLOAT -- maps will be refused")
      else if (span > 1.25 || max(sz) > 1.5 * nf)
        note("membership varies materially across trees")
      else
        note("membership is essentially stable (%d-%d tips) -- the variation is minor",
             min(sz), max(sz))
    }
    if (max(sz) > MAX_FOCAL_BLOAT * nf)
      bad("focal MRCA clade exceeds MAX_FOCAL_BLOAT in at least one tree -- 13_figures.R skips those maps")
  }
}

# 5. families with sequence, and which of them the backbone can constrain
cat("\nfamilies\n")
fam_all <- unique(tx[species %chin% gsub("_", " ", allsp), family])
cat(sprintf("  %d families among the %d sequenced species\n", length(fam_all), length(allsp)))
if (file.exists("backbone_families.tsv")) {
  bb <- unique(fread("backbone_families.tsv")$family)
  nf <- setdiff(fam_all, bb)
  if (length(nf))
    note("%d families are NOT in the backbone, so are unconstrained: %s",
         length(nf), paste(nf, collapse = ", "))
  ok("%d families can be constrained", length(intersect(fam_all, bb)))
}

cat(sprintf("\n%s\n", if (fail) sprintf("%d PROBLEM(S) -- do not submit 08_tree.sbatch yet", fail)
                      else "all checks passed -- safe to submit 08_tree.sbatch"))
quit(status = if (fail) 1L else 0L)
