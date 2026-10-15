#!/usr/bin/env Rscript
# Date the ML tree and the bootstrap trees with LSD2 (least-squares dating, To et al. 2016). Node ages
# from the published backbone chronogram are fixed as secondary calibrations (family-, genus- and
# species-level nodes shared with our tree), and each calibrated clade has its own substitution rate.
# Every dated tree is checked; trees that fail a check are left out and listed in tree/date_failures.tsv.
# usage: Rscript 09_date.R [--ml-only]      (LSD2 binary: on PATH, or LSD2=/path/to/lsd2)

suppressPackageStartupMessages({ library(data.table); library(ape) })
source("config.R")
source("tree_utils.R")
source("root_tree.R")
set.seed(SEED)

LSD2         <- Sys.getenv("LSD2", "lsd2")
MIN_BLEN     <- 0.1    # shortest allowed branch, Ma (LSD2 -u/-U; its default 0 gives zero-length branches)
CAL_LEVELS   <- c("family", "genus", "species")   # backbone nodes used as calibrations (see build_calib)
LOCAL_CLOCKS <- TRUE   # own rate for each calibrated clade (LSD2 -p); FALSE = one rate for the whole tree
ML_ONLY      <- "--ml-only" %in% commandArgs(TRUE)
WORK         <- file.path(tempdir(), "lsd2")
if (!nzchar(Sys.which(LSD2))) stop("lsd2 not found -- put it on PATH or set LSD2=/path/to/lsd2", call. = FALSE)

# alignment length for LSD2's variances (-s), from the IQ-TREE report
iq    <- grep("total sites", readLines("support/support.iqtree"), value = TRUE)[1]
NSITE <- as.integer(sub(".* ([0-9]+) total sites.*", "\\1", iq))
if (is.na(NSITE)) stop("alignment length not found in support/support.iqtree", call. = FALSE)

taxa <- fread(TAXA_TSV); taxa[, tip := gsub(" ", "_", species)]

# backbone chronogram (Newick or NEXUS), node ages in Ma
bb <- if (any(grepl("#NEXUS", readLines(BACKBONE_TRE, n = 5), ignore.case = TRUE)))
        read.nexus(BACKBONE_TRE) else read.tree(BACKBONE_TRE)
bb$tip.label <- gsub("^\'|\'$", "", bb$tip.label)
if (!is.ultrametric(bb, tol = 1e-4)) warning("backbone is not ultrametric")
bb_age <- branching.times(bb) * BB_TIME_SCALE   # see BB_TIME_SCALE in config.R

bbf <- fread("backbone_families.tsv")   # written by 00_resolve_names.R
bbf[, bb_tip := gsub("^'|'$", "", bb_tip)]

# for families and for genera: backbone tip indices per group, the groups below every backbone node,
# and our tips per group
bb_nodes <- (Ntip(bb) + 1):(Ntip(bb) + Nnode(bb))
bb_tips  <- tips_by_node(bb)
grp <- list()
for (lv in setdiff(CAL_LEVELS, "species")) {
  of    <- setNames(bbf[[lv]], bbf$bb_tip)
  below <- list()
  for (nd in bb_nodes) below[[as.character(nd)]] <- unique(na.omit(of[bb$tip.label[bb_tips[[nd]]]]))
  grp[[lv]] <- list(tip_of = split(seq_along(bb$tip.label), of[bb$tip.label]), below = below,
                    ours = split(taxa$tip, taxa[[lv]]))
}

# species-level: backbone tip -> our tip name, for the species in both trees
sp_of <- setNames(gsub(" ", "_", bbf$species), bbf$bb_tip)
sp_of <- sp_of[sp_of %in% taxa$tip]
# genera: of our tips, of backbone tips, and every genus the backbone contains
genus_of_tip <- setNames(taxa$genus, taxa$tip)
bb_genus     <- setNames(bbf$genus, bbf$bb_tip)[bb$tip.label]
bb_genera    <- unique(na.omit(bb_genus))

# calibration table for one tree. Family and genus level: a backbone node's age is used only if that
# node IS the MRCA of the families (genera) it shares with our tree (otherwise its age belongs to a
# deeper split) and those families (genera) are monophyletic in our tree. Species level: every backbone
# node whose species shared with our tree are exactly the shared species below one of our nodes
# (congruent nodes, as in congruification; Eastman et al. 2013) and whose genera agree with that clade's
build_calib <- function(tr) {
  out <- list()
  for (lv in names(grp)) {
    g <- grp[[lv]]
    for (nd in names(g$below)) {
      gs <- intersect(g$below[[nd]], names(g$ours))
      if (length(gs) < 2) next
      if (!identical(as.integer(getMRCA(bb, unlist(g$tip_of[gs]))), as.integer(nd))) next
      tips <- unlist(g$ours[gs], use.names = FALSE)
      tips <- tips[tips %in% tr$tip.label]
      if (length(tips) < 2) next
      if (!is.monophyletic(tr, tips)) next
      a <- unname(bb_age[as.character(nd)])
      if (is.na(a)) next
      out[[length(out) + 1]] <- data.table(node = getMRCA(tr, tips), age.min = a, age.max = a,
                                           soft.bounds = FALSE, level = lv, ngrp = length(gs))
    }
  }
  if ("species" %in% CAL_LEVELS) {
    shared <- intersect(unname(sp_of), tr$tip.label)
    bb_sh  <- names(sp_of)[sp_of %in% shared]
    tr_tips <- tips_by_node(tr)
    for (nd in bb_nodes) {
      b <- intersect(bb$tip.label[bb_tips[[nd]]], bb_sh)
      if (length(b) < 2) next
      if (getMRCA(bb, b) != nd) next
      s <- unique(unname(sp_of[b]))
      if (length(s) < 2) next
      m <- getMRCA(tr, s)
      if (!setequal(intersect(tr$tip.label[tr_tips[[m]]], shared), s)) next
      # the genera must agree too: no genus the backbone places outside this node may sit inside our
      # clade, and no genus inside the backbone node may have members of ours outside the clade
      in_ours <- unique(genus_of_tip[tr$tip.label[tr_tips[[m]]]])
      in_bb   <- unique(na.omit(bb_genus[bb_tips[[nd]]]))
      if (length(setdiff(intersect(in_ours, bb_genera), in_bb))) next
      out_ours <- unique(genus_of_tip[tr$tip.label[-tr_tips[[m]]]])
      if (length(intersect(intersect(in_bb, out_ours), taxa$genus))) next
      a <- unname(bb_age[as.character(nd)])
      if (is.na(a)) next
      out[[length(out) + 1]] <- data.table(node = m, age.min = a, age.max = a,
                                           soft.bounds = FALSE, level = "species", ngrp = length(s))
    }
  }
  if (!length(out)) stop("no calibratable nodes -- check family name matching")
  cal <- rbindlist(out)
  # one calibration per node: family before genus before species, then the most inclusive
  cal <- cal[order(node, match(level, c("family", "genus", "species")), -ngrp)][!duplicated(node)]

  # drop genus- and species-level calibrations that are not younger than every calibrated ancestor
  # and older than every calibrated descendant (the trees can disagree)
  anc  <- ancestors_by_node(tr)
  drop <- logical(nrow(cal))
  for (i in which(cal$level != "family")) for (j in seq_len(nrow(cal))) {
    if (i == j) next
    if (cal$node[j] %in% anc[[cal$node[i]]] && cal$age.min[j] <= cal$age.min[i]) drop[i] <- TRUE
    if (cal$node[i] %in% anc[[cal$node[j]]] && cal$age.min[j] >= cal$age.min[i]) drop[i] <- TRUE
  }
  attr_n <- c(family = sum(cal$level == "family"), genus = sum(cal$level == "genus" & !drop),
              species = sum(cal$level == "species" & !drop), dropped = sum(drop))
  cal <- cal[!drop]
  cal[, ngrp := NULL]
  cal <- as.data.frame(cal)
  attr(cal, "n") <- attr_n
  cal
}

# root on the outgroup suborder, resolve, and drop node labels (support values) before dating
prep_tree <- function(tr) {
  tr$tip.label <- sub("^_R_", "", tr$tip.label)   # mafft --adjustdirection prefix
  if (!is.rooted(tr)) tr <- root_tree(tr, taxa)
  tr <- multi2di(tr)
  tr$edge.length[tr$edge.length <= 0] <- 1e-8
  tr$node.label <- NULL
  tr
}

# reason a dated tree is unusable, or NULL: tips, rooted topology, branch lengths, ultrametric, calibrations
check_dated <- function(dt, tr, cal, pair) {
  dt$node.label <- NULL
  if (!setequal(dt$tip.label, tr$tip.label)) return("tip set changed")
  if (!is.rooted(dt) || Nnode(dt) != Ntip(dt) - 1) return("not rooted and binary")
  if (!isTRUE(all.equal.phylo(dt, tr, use.edge.length = FALSE))) return("topology changed")
  bad <- dt$edge.length <= 0
  if (any(bad)) return(sprintf("%d branches of length <= 0 (%d to tips, %d exactly 0, shortest %.3g)",
                               sum(bad), sum(bad & dt$edge[, 2] <= Ntip(dt)),
                               sum(dt$edge.length == 0), min(dt$edge.length)))
  dep <- node.depth.edgelength(dt)[seq_len(Ntip(dt))]
  if (diff(range(dep)) > 1e-4 * max(dep)) return("not ultrametric")
  bt <- unname(branching.times(dt))
  for (k in seq_len(nrow(cal))) {
    a <- bt[getMRCA(dt, pair[k, ]) - Ntip(dt)]
    if (abs(a - cal$age.min[k]) > 1e-3 * cal$age.min[k])
      return(sprintf("calibrated node %s,%s at %.2f Ma, not %.2f", pair[k, 1], pair[k, 2], a, cal$age.min[k]))
  }
  NULL
}

# date one tree in its own folder; returns the dated tree (NULL if rejected), the reason and LSD2's output.
# min_blen: shortest allowed branch in Ma (LSD2 -u and -U); LSD2's default is 0
date_tree <- function(tr, dir, min_blen = 0, cal = NULL) {
  if (is.null(cal)) { tr <- prep_tree(tr); cal <- build_calib(tr) }

  # each calibrated node named as the MRCA of one tip from each of its two daughter clades
  tips <- tips_by_node(tr)
  pair <- matrix("", nrow(cal), 2)
  for (k in seq_len(nrow(cal))) {
    kids <- tr$edge[tr$edge[, 1] == cal$node[k], 2]
    pair[k, ] <- tr$tip.label[c(tips[[kids[1]]][1], tips[[kids[2]]][1])]
  }

  # tree and date file; ages in Ma become negative dates, every tip is dated 0 (-z 0)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  tf <- file.path(dir, "in.nwk"); df <- file.path(dir, "dates.txt"); out <- file.path(dir, "out")
  writeLines(c(nrow(cal), sprintf("mrca(%s,%s) %.6f", pair[, 1], pair[, 2], -cal$age.min)), df)

  # local clocks: every calibrated clade gets its own rate, from its root down to the next calibrated
  # clades below it; branches outside every calibrated clade share one more rate
  part <- character(0)
  if (LOCAL_CLOCKS) {
    tl <- tr; tl$node.label <- rep("", Nnode(tr))
    tl$node.label[cal$node - Ntip(tr)] <- sprintf("cal%d", seq_len(nrow(cal)))
    anc <- ancestors_by_node(tr)
    lines <- character(0)
    for (k in seq_len(nrow(cal))) {
      if (cal$node[k] == Ntip(tr) + 1) next   # the root's clade is the remaining-branches group
      inner <- integer(0)
      for (j in seq_len(nrow(cal))) {
        if (j == k || !(cal$node[k] %in% anc[[cal$node[j]]])) next
        between <- intersect(anc[[cal$node[j]]], cal$node)
        if (cal$node[k] == between[1]) inner <- c(inner, j)   # nearest calibrated ancestor of j is k
      }
      lines <- c(lines, sprintf("r%d 1 {%s}", k, paste(sprintf("cal%d", c(k, inner)), collapse = " ")))
    }
    part <- c("-p", file.path(dir, "rates.txt"))
    write.tree(tl, tf); writeLines(lines, part[2])
  } else write.tree(tr, tf)

  # LSD2: keep the root on the outgroup branch but place it along that branch (-r k),
  # collapse no short branches (-l -1)
  log <- suppressWarnings(system2(LSD2, c("-i", tf, "-d", df, "-s", NSITE, "-z", 0, "-l", -1,
                                          "-u", min_blen, "-U", min_blen, part, "-r", "k", "-o", out),
                                  stdout = TRUE, stderr = TRUE))
  res <- paste0(out, ".date.nexus")
  if (!file.exists(res))
    return(list(tree = NULL, why = paste("no", basename(res), "from LSD2:", tail(c("", log), 1)), log = log))

  # time-scaled tree, with LSD2's [&...] node annotations removed
  nx <- gsub("\\[[^]]*\\]", "", readLines(res))
  writeLines(nx, paste0(out, ".plain.nexus"))
  dt <- read.nexus(paste0(out, ".plain.nexus"))
  if (inherits(dt, "multiPhylo")) dt <- dt[[1]]
  dt$node.label <- NULL
  dt$tip.label <- gsub("^'|'$", "", dt$tip.label)

  why <- check_dated(dt, tr, cal, pair)
  list(tree = if (is.null(why)) dt else NULL, why = why, log = log, ncal = attr(cal, "n"))
}

# ML tree; its LSD2 input and output files are kept in tree/lsd2_ml/ (LSD2's report is tree/lsd2_ml/out)
cat(sprintf("LSD2: %s; alignment length %d; minimum branch length %g Ma\n", Sys.which(LSD2), NSITE, MIN_BLEN))
ml <- read.tree("support/support.treefile")
r  <- date_tree(ml, "tree/lsd2_ml", MIN_BLEN)
writeLines(r$log, "tree/lsd2_ml/lsd2_screen.log")
cat(sprintf("calibrations on the ML tree: %d family-level, %d genus-level, %d species-level (%d dropped as inconsistent)%s\n",
            r$ncal["family"], r$ncal["genus"], r$ncal["species"], r$ncal["dropped"],
            if (LOCAL_CLOCKS) "; own rate per calibrated clade" else ""))
if (is.null(r$tree)) stop("ML tree rejected: ", r$why, "\n  (LSD2 files in tree/lsd2_ml/)", call. = FALSE)
ml_dated <- r$tree
write.tree(ml_dated, "tree/odon.dated.tre")

# summary of the dated ML tree
foc <- intersect(gsub(" ", "_", FOCAL_TIPS), ml_dated$tip.label)
bt  <- unname(branching.times(ml_dated))
cat(sprintf("ML tree dated: root %.1f Ma, %s crown %.2f Ma, shortest branch %.4f Ma (%d branches at the minimum)\n",
            max(bt), FOCAL_GENUS, bt[getMRCA(ml_dated, foc) - Ntip(ml_dated)], min(ml_dated$edge.length),
            sum(ml_dated$edge.length < MIN_BLEN * 1.001)))
if (ML_ONLY) quit(save = "no")

# bootstrap trees
bs <- read.tree("support/support.ufboot")
if (inherits(bs, "phylo")) bs <- list(bs)
if (!length(bs)) stop("support.ufboot is empty -- 08b_support.sbatch needs --wbtl")
cat(sprintf("dating %d bootstrap trees on %d cores\n", length(bs), MC_CORES))
t0 <- Sys.time()
res <- parallel::mclapply(seq_along(bs), function(k) {
  d <- file.path(WORK, k)
  r <- try(date_tree(bs[[k]], d, MIN_BLEN), silent = TRUE)
  unlink(d, recursive = TRUE)
  if (inherits(r, "try-error")) list(tree = NULL, why = conditionMessage(attr(r, "condition"))) else r[c("tree", "why")]
}, mc.cores = MC_CORES)
cat(sprintf("  done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# keep the trees that passed every check; list the rest with the reason
ok    <- !vapply(res, function(r) is.null(r$tree), logical(1))
dated <- lapply(res[ok], `[[`, "tree")
class(dated) <- "multiPhylo"
write.tree(dated, "tree/odon.boot.dated.tre")
fails <- data.table(tree = which(!ok), reason = vapply(res[!ok], function(r) as.character(r$why), ""))
fwrite(fails, "tree/date_failures.tsv", sep = "\t")

# root and focal crown ages over the dated bootstrap trees
root_age <- numeric(0); crown_age <- numeric(0)
for (d in dated) {
  b <- unname(branching.times(d))
  root_age  <- c(root_age, max(b))
  crown_age <- c(crown_age, b[getMRCA(d, intersect(foc, d$tip.label)) - Ntip(d)])
}
cat(sprintf("%d dated, %d rejected (tree/date_failures.tsv)\n", sum(ok), sum(!ok)))
if (sum(ok)) cat(sprintf("root age %.1f-%.1f Ma; %s crown median %.2f Ma (%.2f-%.2f)\n",
                         min(root_age), max(root_age), FOCAL_GENUS, median(crown_age),
                         min(crown_age), max(crown_age)))
if (nrow(fails)) print(fails[, .N, by = .(reason = sub(" \\(.*", "", reason))])
if (mean(!ok) > 0.05) warning(sprintf("%.1f%% of bootstrap trees rejected -- check tree/date_failures.tsv",
                                      100 * mean(!ok)), call. = FALSE)
