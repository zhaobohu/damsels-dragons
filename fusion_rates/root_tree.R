# Root an (unrooted) IQ-TREE tree on the OUTGROUP_SUBORDER tips of taxa.tsv
suppressPackageStartupMessages({ library(ape); library(data.table) })

root_tree <- function(tr, taxa = NULL, quiet = TRUE) {
  tr$tip.label <- sub("^_R_", "", tr$tip.label)
  if (is.null(taxa)) taxa <- fread(TAXA_TSV)
  og <- gsub(" ", "_", taxa[suborder == OUTGROUP_SUBORDER, species])
  og <- intersect(og, tr$tip.label)

  # fewer than 2 outgroup tips: leave unrooted, with a warning
  if (length(og) < 2) {
    warning("outgroup ", OUTGROUP_SUBORDER, ": ", length(og),
            " tips present -- tree left unrooted"); return(tr)
  }

  # monophyletic outgroup: root on all of it
  if (is.monophyletic(tr, og)) {
    rt <- root(tr, outgroup = og, resolve.root = TRUE)
    if (!quiet) message(sprintf("rooted on %d %s tips", length(og), OUTGROUP_SUBORDER))
    return(rt)
  }

  # otherwise root on a single outgroup tip (any root outside the focal clade will do)
  rt <- root(tr, outgroup = og[1], resolve.root = TRUE)
  warning("outgroup ", OUTGROUP_SUBORDER, " is not monophyletic (",
          length(og), " tips) -- rooted on a single tip, ", og[1],
          ". Check the tree before trusting clade ages.")
  rt
}
