# Base-R replacements for phangorn::Descendants and phangorn::Ancestors (avoids the fastmatch dependency)

suppressPackageStartupMessages(library(ape))

# tip indices below every node, tips included (a tip subtends itself)
tips_by_node <- function(tr) {
  ntip <- Ntip(tr); nn <- ntip + tr$Nnode
  out <- vector("list", nn)
  for (i in seq_len(ntip)) out[[i]] <- i
  # postorder: every child is finished before its parent
  e <- reorder(tr, "postorder")$edge
  for (k in seq_len(nrow(e))) out[[e[k, 1L]]] <- c(out[[e[k, 1L]]], out[[e[k, 2L]]])
  lapply(out, function(v) if (is.null(v)) integer(0) else sort.int(v))
}

# ancestors of every node, from parent up to the root (root gets integer(0))
ancestors_by_node <- function(tr) {
  nn <- Ntip(tr) + tr$Nnode
  out <- vector("list", nn)
  for (i in seq_len(nn)) out[[i]] <- integer(0)
  # cladewise: every parent comes before its children
  e <- reorder(tr, "cladewise")$edge
  for (k in seq_len(nrow(e))) {
    pa <- e[k, 1L]; ch <- e[k, 2L]
    out[[ch]] <- c(pa, out[[pa]])
  }
  out
}
