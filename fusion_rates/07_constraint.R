#!/usr/bin/env Rscript
# Family-level IQ-TREE topology constraint from the dated backbone, and the partition file

suppressPackageStartupMessages({ library(data.table); library(ape) })
source("config.R")

# taxa minus drop_tips.txt, restricted to species with sequence in any marker (tree_taxa.txt)
taxa <- fread(TAXA_TSV)
drop <- if (file.exists("drop_tips.txt")) readLines("drop_tips.txt") else character(0)
drop <- drop[nzchar(drop)]
taxa <- taxa[!gsub(" ", "_", species) %in% drop]
taxa[, tip := gsub(" ", "_", species)]
taxa <- taxa[tip %in% sub("^_R_", "", readLines("tree_taxa.txt"))]

# backbone (nexus or newick), tip quotes stripped
bb <- if (any(grepl("#NEXUS", readLines(BACKBONE_TRE, n = 5), ignore.case = TRUE)))
        read.nexus(BACKBONE_TRE) else read.tree(BACKBONE_TRE)
bb$tip.label <- gsub("^\'|\'$", "", bb$tip.label)

# backbone tip -> family
bbf <- fread("backbone_families.tsv")   # written by 00_resolve_names.R
bbf[, bb_tip := gsub("^'|'$", "", bb_tip)]

# families shared with the backbone; the rest are left unconstrained
fams <- intersect(unique(taxa$family), unique(bbf$family))
missing <- setdiff(unique(taxa$family), fams)
if (length(missing))
  warning("families absent from backbone, left unconstrained: ",
          paste(missing, collapse = ", "))

# one backbone exemplar per family gives the inter-family topology
ex <- bbf[, .(bb_tip = bb_tip[1]), by = family][family %in% fams]
sk <- keep.tip(bb, ex$bb_tip)
sk$tip.label <- ex$family[match(sk$tip.label, ex$bb_tip)]

# replace each family tip with a polytomy of its species
nw <- write.tree(sk)
for (f in fams) {
  sp <- taxa[family == f, tip]
  rep <- if (length(sp) == 1) sp else sprintf("(%s)", paste(sp, collapse = ","))
  nw <- sub(sprintf("(?<![A-Za-z_])%s(?![A-Za-z_])", f), rep, nw, perl = TRUE)
}
writeLines(nw, "tree/constraint.tre")

# partition file: one charset per codon position for coding markers, one otherwise
coding <- c("COI", "EF1a", "H3")
lines <- c("#nexus", "begin sets;")
k <- 0
present <- sub("\\.trim\\.fasta$", "", basename(
  list.files("align", "\\.trim\\.fasta$")))
for (m in intersect(MARKERS, present)) {
  f <- sprintf("align/%s.trim.fasta", m)
  if (!file.exists(f) || file.size(f) == 0) next
  a <- readLines(f); h <- grep("^>", a)
  L <- sum(nchar(a[(h[1] + 1):(if (length(h) > 1) h[2] - 1 else length(a))]))
  if (m %in% coding) {
    for (cp in 1:3) {
      k <- k + 1
      lines <- c(lines, sprintf("  charset p%d = %s: %d-%d\\3;", k, f, cp, L))
    }
  } else {
    k <- k + 1
    lines <- c(lines, sprintf("  charset p%d = %s: 1-%d;", k, f, L))
  }
}
lines <- c(lines, "end;")
writeLines(lines, "tree/partitions.nex")

cat(sprintf("constraint: %d families, %d tips\n", length(fams), nrow(taxa)))
cat(sprintf("partitions: %d\n", k))
