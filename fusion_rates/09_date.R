#!/usr/bin/env Rscript
# Date the ML tree and the bootstrap trees with chronos, fixing node ages read off the published
# backbone chronogram as secondary calibrations (min = max)

suppressPackageStartupMessages({ library(data.table); library(ape) })
source("config.R")
source("tree_utils.R")
source("root_tree.R")
set.seed(SEED)

taxa <- fread(TAXA_TSV); taxa[, tip := gsub(" ", "_", species)]

# backbone chronogram (Newick or NEXUS), node ages in Ma
bb <- if (any(grepl("#NEXUS", readLines(BACKBONE_TRE, n = 5), ignore.case = TRUE)))
        read.nexus(BACKBONE_TRE) else read.tree(BACKBONE_TRE)
bb$tip.label <- gsub("^\'|\'$", "", bb$tip.label)
if (!is.ultrametric(bb, tol = 1e-4)) warning("backbone is not ultrametric")
bb_age <- branching.times(bb) * BB_TIME_SCALE   # see BB_TIME_SCALE in config.R

bbf <- fread("backbone_families.tsv")   # written by 00_resolve_names.R
bbf[, bb_tip := gsub("^'|'$", "", bb_tip)]

# families below every backbone node, and backbone tip indices per family
bb_nodes <- (Ntip(bb) + 1):(Ntip(bb) + Nnode(bb))
fam_of   <- setNames(bbf$family, bbf$bb_tip)
bb_tips  <- tips_by_node(bb)
bb_tip_of <- split(seq_along(bb$tip.label), fam_of[bb$tip.label])
bb_fams <- list()
for (nd in bb_nodes) {
  tips <- bb$tip.label[bb_tips[[nd]]]
  bb_fams[[as.character(nd)]] <- unique(na.omit(fam_of[tips]))
}

# calibration table for one tree: a backbone node's age is used only if that node IS the MRCA
# of the families it shares with our tree (otherwise its age belongs to a deeper split)
# and those families are monophyletic in our tree
build_calib <- function(tr) {
  fam_tip <- split(taxa$tip, taxa$family)
  out <- list()
  for (nd in names(bb_fams)) {
    fs <- intersect(bb_fams[[nd]], names(fam_tip))
    if (length(fs) < 2) next
    if (!identical(as.integer(getMRCA(bb, unlist(bb_tip_of[fs]))), as.integer(nd))) next
    tips <- unlist(fam_tip[fs], use.names = FALSE)
    tips <- tips[tips %in% tr$tip.label]
    if (length(tips) < 2) next
    if (!is.monophyletic(tr, tips)) next
    mr <- getMRCA(tr, tips)
    a  <- unname(bb_age[as.character(nd)])
    if (is.na(a)) next
    out[[length(out) + 1]] <- data.table(node = mr, age.min = a, age.max = a,
                                         soft.bounds = FALSE, nfam = length(fs))
  }
  if (!length(out)) stop("no calibratable nodes -- check family name matching")
  cal <- rbindlist(out)
  # keep the deepest (most inclusive) calibration per node
  cal <- cal[order(node, -nfam)][!duplicated(node)]
  cal[, nfam := NULL]
  as.data.frame(cal)
}

# root (chronos needs a rooted tree), resolve, calibrate and date one tree; NULL if chronos fails
date_tree <- function(tr, lambda = 1) {
  if (!is.rooted(tr)) tr <- root_tree(tr, taxa)
  tr <- multi2di(tr)
  tr$edge.length[tr$edge.length <= 0] <- 1e-8
  cal <- build_calib(tr)
  ct <- try(chronos(tr, lambda = lambda, model = CHRONOS_MODEL,
                    calibration = cal, quiet = TRUE), silent = TRUE)
  if (inherits(ct, "try-error")) return(NULL)
  class(ct) <- "phylo"; ct
}

# ML tree dated at each lambda; keep the fit with the lowest PHIIC
ml <- read.tree("support/support.treefile")
ml$tip.label <- sub("^_R_", "", ml$tip.label)   # mafft --adjustdirection* prefix
lams <- c(0.1, 1, 10, 100)
fits <- vector("list", length(lams))
for (k in seq_along(lams)) fits[k] <- list(date_tree(ml, lams[k]))
ok   <- !vapply(fits, is.null, logical(1))
phi  <- numeric(0)
for (f in fits[ok]) phi <- c(phi, attr(f, "PHIIC")$PHIIC)
best <- lams[ok][which.min(phi)]
cat("lambda by PHIIC:\n"); print(data.frame(lambda = lams[ok], PHIIC = phi))
cat(sprintf("chosen lambda = %g\n", best))

ml_dated <- fits[ok][[which.min(phi)]]
write.tree(ml_dated, "tree/odon.dated.tre")

# bootstrap trees
bs <- read.tree("support/support.ufboot")
if (!length(bs)) stop("support.ufboot is empty -- 08b_support.sbatch needs --wbtl")
for (k in seq_along(bs)) bs[[k]]$tip.label <- sub("^_R_", "", bs[[k]]$tip.label)

set.seed(SEED)
idx <- seq_along(bs)
cat(sprintf("dating %d bootstrap trees\n", length(idx)))

# resume from the checkpoint if there is one
CKPT <- "tree/boot_dated_ckpt.rds"
done <- if (file.exists(CKPT)) readRDS(CKPT) else list()
if (length(done)) cat(sprintf("  resuming: %d already dated in %s\n", length(done), CKPT))
todo <- setdiff(as.character(idx), names(done))

# date the rest in parallel chunks, checkpointing after each chunk
chunk <- max(1L, MC_CORES * 4L)
t0 <- Sys.time()
starts <- if (length(todo)) seq(1, length(todo), by = chunk) else integer(0)
if (!length(starts)) cat("  nothing left to date\n")
for (start in starts) {
  part <- todo[start:min(start + chunk - 1L, length(todo))]
  res  <- parallel::mclapply(part, function(k) date_tree(bs[[as.integer(k)]], best),
                             mc.cores = MC_CORES)
  names(res) <- part
  done <- c(done, res)
  saveRDS(done, CKPT)
  el  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  n   <- start + length(part) - 1L
  eta <- el / n * (length(todo) - n)
  cat(sprintf("  %d/%d dated  (%.0f min elapsed, ~%.0f min left)\n",
              n, length(todo), el / 60, eta / 60))
  flush.console()
}

# write the dated bootstrap trees (failures dropped) and report root ages
dated <- done[as.character(idx)]
fail  <- vapply(dated, is.null, logical(1))
dated <- dated[!fail]
class(dated) <- "multiPhylo"
write.tree(dated, "tree/odon.boot.dated.tre")
cat(sprintf("%d dated (%d failed)\n", length(dated), sum(fail)))
cat(sprintf("checkpoint kept at %s -- delete it to redate from scratch\n", CKPT))
root_age <- numeric(0)
for (d in dated) root_age <- c(root_age, max(branching.times(d)))
cat(sprintf("root age: ML %.1f Ma, bootstrap %.1f-%.1f Ma\n",
            max(branching.times(ml_dated)), min(root_age), max(root_age)))
