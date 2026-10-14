library(data.table)
library(viridis)

# anchor-partner rank filter, shared with divergence_dating/02.select_loci.R and figs/fig4/fig4.R
source("../../lib/anchor_filter.R")

# the pairs.sbatch run directory (gene_pairs_ds)
run_dir <- "../../gene_pairs_ds"

# chromosomes (bp = hard-coded fallback fusion breakpoint on the Aeshna X)
ENT <- list(
  juncea_X = list(name = "Aeshna X", len = 184464265,
                  bp = 103.4e6, col_neo = "#7570B3", col_anc = "#66A61E"),
  juncea_Y = list(name = "Aeshna Y", len = 166739284, col = "#B3717C"),
  anax_12  = list(name = "Anax 12",  len = 76569406,  col = "#7570B3"),
  anax_X   = list(name = "Anax X",   len = 73439698,  col = "#66A61E")
)

# quality filters (min.ratio / min.pid / ratio.na.keep are also read by anchor_filter.R)
min.score <- 0.5
min.frac  <- 0.3
cov.def   <- 0.5    # min(query, subject) coverage to call a union pair 1:1; 0 = raw union
min.ratio <- 0.65   # anchor partner rank floor; 0 = no rank filter
min.pid   <- 0      # identity floor ANDed with min.ratio; 0 = rank alone
ratio.na.keep <- FALSE  # drop anchor pairs with no row in the rank search

# margins + how far the side annotations sit from the figure
fig.mar  <- c(3, 12, 3, 12)
lab.frac <- 0.02

strip_tx <- function(x) sub("\\.t[0-9]+$", "", x)

# gene midpoints from a chromosome's longest-isoform GFF
readGenes <- function(label) {
  path <- file.path(run_dir, "chrom", paste0(label, ".longest.gff3"))
  if (!file.exists(path)) { warning("missing gff: ", path); return(NULL) }
  gff <- read.delim(path, comment.char = "#", header = FALSE)
  gff <- gff[gff$V3 == "gene", ]
  gff$gene_id <- gsub(".*ID=([^;]+).*", "\\1", gff$V9)
  gff$mid     <- rowMeans(gff[, 4:5])
  gff[, c("gene_id", "mid")]
}

bezierCubic <- function(p0, p1, p2, p3, n = 100) {
  t <- seq(0, 1, length.out = n)
  list(
    x = (1-t)^3*p0[1] + 3*(1-t)^2*t*p1[1] + 3*(1-t)*t^2*p2[1] + t^3*p3[1],
    y = (1-t)^3*p0[2] + 3*(1-t)^2*t*p1[2] + 3*(1-t)*t^2*p2[2] + t^3*p3[2]
  )
}

# union pairs of one comparison (comp/<cmp>/pairs.union.tsv), optionally only those with
# coverage >= cov_min
readUnion <- function(cmp, cov_min = 0) {
  p <- file.path(run_dir, "comp", cmp, "pairs.union.tsv")
  if (!file.exists(p)) { warning("missing union: ", p); return(NULL) }
  u <- as.data.frame(fread(p, header = FALSE))[, 1:2]
  names(u) <- c("gene_a", "gene_b")
  if (cov_min > 0 && !is.null(covtab)) {
    cv <- covtab[covtab$comparison == cmp, c("gene_a", "gene_b", "min_cov")]
    u  <- merge(u, cv, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
    u  <- u[!is.na(u$min_cov) & u$min_cov >= cov_min, c("gene_a", "gene_b"), drop = FALSE]
  }
  u
}

# 1:1 pairs of one comparison: reciprocal singletons of the union, with dS/quality columns
# attached where present
one2one <- function(cmp, tbl = dat, cov_min = cov.def) {
  u <- readUnion(cmp, cov_min = cov_min)
  if (is.null(u) || !nrow(u)) return(tbl[tbl$comparison == cmp, , drop = FALSE][0, ])
  oo <- u[ave(rep(1, nrow(u)), u$gene_a, FUN = length) == 1 &
            ave(rep(1, nrow(u)), u$gene_b, FUN = length) == 1, ]
  q  <- tbl[tbl$comparison == cmp, setdiff(names(tbl), "comparison"), drop = FALSE]
  m  <- merge(oo, q, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
  m$comparison <- cmp
  m
}

# attach the anchor partner rank (anchor_filter.R) and drop anchors below the floors
add_ratio <- function(d, cmp, ratio_min = min.ratio, pid_min = min.pid,
                      na_keep = ratio.na.keep) {
  attach_rank(d, cmp, run_dir, ratio_min = ratio_min,
              pid_min = pid_min, na_keep = na_keep)
}

# MCScanX anchor partner of each top-chromosome gene (comp/<cmp>/mcscanx.pairs.tsv);
# one partner per top gene: best coverage, then lowest dS
synteny_partner <- function(top, cmps, ratio_min = min.ratio) {
  rows <- do.call(rbind, lapply(cmps, function(cmp) {
    parts <- strsplit(cmp, "__vs__")[[1]]
    p   <- file.path(run_dir, "comp", cmp, "mcscanx.pairs.tsv")
    ida <- file.path(run_dir, "chrom", paste0(parts[1], ".ids"))
    if (!file.exists(p) || file.size(p) == 0 || !file.exists(ida)) {
      warning("no MCScanX anchors for ", cmp, " (", p, ")"); return(NULL)
    }
    m    <- as.data.frame(fread(p, header = FALSE))[, 1:2]; names(m) <- c("g1", "g2")
    aids <- readLines(ida)
    a1   <- m$g1 %in% aids                                   # orient: gene_a = parts[1] side
    d <- data.frame(comparison = cmp,
                    gene_a = ifelse(a1, m$g1, m$g2),
                    gene_b = ifelse(a1, m$g2, m$g1), stringsAsFactors = FALSE)
    d$top.gene <- strip_tx(if (top == parts[1]) d$gene_a else d$gene_b)
    if (!is.null(covtab)) {
      cv <- covtab[covtab$comparison == cmp, c("gene_a", "gene_b", "min_cov")]
      d  <- merge(d, cv, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
    } else d$min_cov <- NA_real_
    ds <- dat[dat$comparison == cmp, c("gene_a", "gene_b", "dS")]
    d  <- merge(d, ds, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
    add_ratio(d, cmp, ratio_min)
  }))
  if (is.null(rows) || !nrow(rows)) return(rows)
  s  <- ifelse(is.na(rows$min_cov), 0, rows$min_cov)
  dd <- ifelse(is.na(rows$dS), Inf, rows$dS)
  rows <- rows[order(rows$top.gene, -s, dd), ]
  rows[!duplicated(rows$top.gene), , drop = FALSE]
}

# MCScanX anchors that are reciprocal singletons (1:1 bijection)
synteny_recip_partner <- function(top, cmps, ratio_min = min.ratio) {
  rows <- do.call(rbind, lapply(cmps, function(cmp) {
    parts <- strsplit(cmp, "__vs__")[[1]]
    p   <- file.path(run_dir, "comp", cmp, "mcscanx.pairs.tsv")
    ida <- file.path(run_dir, "chrom", paste0(parts[1], ".ids"))
    if (!file.exists(p) || file.size(p) == 0 || !file.exists(ida)) {
      warning("no MCScanX anchors for ", cmp, " (", p, ")"); return(NULL)
    }
    m    <- as.data.frame(fread(p, header = FALSE))[, 1:2]; names(m) <- c("g1", "g2")
    aids <- readLines(ida)
    a1   <- m$g1 %in% aids                                  # g1 sits on parts[1] side
    d <- unique(data.frame(gene_a = ifelse(a1, m$g1, m$g2),
                           gene_b = ifelse(a1, m$g2, m$g1),
                           stringsAsFactors = FALSE))
    ka <- strip_tx(d$gene_a); kb <- strip_tx(d$gene_b)
    keep <- ave(rep(1, nrow(d)), ka, FUN = length) == 1 &
            ave(rep(1, nrow(d)), kb, FUN = length) == 1     # reciprocal singleton -> bijection
    d <- d[keep, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d$comparison <- cmp
    d$top.gene   <- strip_tx(if (top == parts[1]) d$gene_a else d$gene_b)
    if (!is.null(covtab)) {
      cv <- covtab[covtab$comparison == cmp, c("gene_a", "gene_b", "min_cov")]
      d  <- merge(d, cv, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
    } else d$min_cov <- NA_real_
    ds <- dat[dat$comparison == cmp, c("gene_a", "gene_b", "dS")]
    d  <- merge(d, ds, by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)
    add_ratio(d, cmp, ratio_min)
  }))
  rows
}

# draw a chromosome bar at vertical centre yc; split by breakpoint if present
draw_bar <- function(E, yc, off_mb = 0, half_h = 0.06) {
  len_mb <- E$len / 1e6
  if (!is.null(E$bp)) {
    bp_mb <- E$bp / 1e6
    rect(off_mb,          yc - half_h, off_mb + bp_mb,  yc + half_h,
         col = E$col_neo, border = "black", lwd = 2)
    rect(off_mb + bp_mb,  yc - half_h, off_mb + len_mb, yc + half_h,
         col = E$col_anc, border = "black", lwd = 2)
  } else {
    rect(off_mb, yc - half_h, off_mb + len_mb, yc + half_h,
         col = E$col, border = "black", lwd = 2)
  }
}

# all pairs with dN/dS
dat <- as.data.frame(fread(file.path(run_dir, "all_ds_results.tsv")))

# per-pair coverage (comparison, gene_a, gene_b, cov_a, cov_b, min_cov)
cov_path <- file.path(run_dir, "all_pairs_cover.tsv")
covtab <- if (file.exists(cov_path)) as.data.frame(fread(cov_path)) else {
  warning("missing ", cov_path, " -- coverage filter disabled (1:1 falls back to raw union)"); NULL
}

# GENESPACE pangenes (orthogroup pairing, used only for the count comparison below)
PG_PATH <- "../../odonate_synteny/genespace/pangenes/Aeshna_juncea_pangenes.txt.gz"
OG_CONF <- c("PASS", "array")
PG_MAP <- list(
  juncea_X = list(genome = "Aeshna_juncea",        chr = "OZ369554.1"),
  juncea_Y = list(genome = "Aeshna_juncea",        chr = "OZ369555.1"),
  anax_12  = list(genome = "Anax_junius",          chr = "OZ349811.1"),
  anax_X   = list(genome = "Anax_junius",          chr = "OZ349813.1")
)

pg_conf <- NULL
if (file.exists(PG_PATH)) {
  pg <- fread(PG_PATH)
  pg[, g := strip_tx(id)]
  pg[, start := as.numeric(start)]
  if ("pgID" %in% names(pg)) pg[is.na(og), og := -as.numeric(pgID)] else pg[is.na(og), og := -.I]
  pg_conf <- if (is.null(OG_CONF)) pg else pg[flag %chin% OG_CONF]
} else {
  warning("missing pangenes: ", PG_PATH, " -- orthogroup vs MCScanX comparison skipped")
}

# orthogroup pairs of one comparison (reciprocal 1:1 if recip), with dS attached
og_partner_one <- function(cmp, top, recip = TRUE) {
  if (is.null(pg_conf)) return(NULL)
  parts <- strsplit(cmp, "__vs__")[[1]]
  ta <- PG_MAP[[parts[1]]]; tb <- PG_MAP[[parts[2]]]
  if (is.null(ta) || is.null(tb)) { warning("no pangenes map for ", cmp); return(NULL) }
  A <- pg_conf[genome == ta$genome & chr == ta$chr, .(og, idA = id, gA = g, posA = start)]
  B <- pg_conf[genome == tb$genome & chr == tb$chr, .(og, idB = id, gB = g, posB = start)]
  if (!nrow(A) || !nrow(B)) { warning("no pangenes genes for ", cmp); return(NULL) }
  P <- merge(A, B, by = "og", allow.cartesian = TRUE)
  if (!nrow(P)) return(NULL)
  if (recip) {
    P[, `:=`(nA = uniqueN(gA), nB = uniqueN(gB)), by = og]
    P <- P[nA == 1L & nB == 1L]
  }
  if (!nrow(P)) return(NULL)
  top_is_a <- (top == parts[1])
  d <- data.frame(
    comparison = cmp,
    gene_a   = P$idA,  gene_b = P$idB,
    top.gene = if (top_is_a) P$gA   else P$gB,
    top.mid  = if (top_is_a) P$posA else P$posB,
    bot.mid  = if (top_is_a) P$posB else P$posA,
    stringsAsFactors = FALSE)
  ds <- dat[dat$comparison == cmp, c("gene_a", "gene_b", "dS")]
  if (nrow(ds)) {
    ds$ka <- strip_tx(ds$gene_a); ds$kb <- strip_tx(ds$gene_b)
    ds <- ds[!duplicated(ds[c("ka", "kb")]), c("ka", "kb", "dS")]
    d$ka <- strip_tx(d$gene_a);  d$kb <- strip_tx(d$gene_b)
    d <- merge(d, ds, by = c("ka", "kb"), all.x = TRUE, sort = FALSE)
    d$ka <- NULL; d$kb <- NULL
  } else d$dS <- NA_real_
  d
}

# print reciprocal 1:1 pair counts, orthogroup vs MCScanX anchors
if (!is.null(pg_conf)) {
  cat("\n[pairing comparison] reciprocal 1:1 pair counts (orthogroup vs MCScanX anchors):\n")
  for (cmp in c("juncea_X__vs__anax_12", "juncea_X__vs__anax_X", "juncea_X__vs__juncea_Y")) {
    t1 <- strsplit(cmp, "__vs__")[[1]][1]
    no <- og_partner_one(cmp, t1, recip = TRUE);      no <- if (is.null(no)) 0L else nrow(no)
    nm <- synteny_recip_partner(t1, cmp);             nm <- if (is.null(nm)) 0L else nrow(nm)
    cat(sprintf("  %-28s orthogroup %4d   MCScanX %4d\n", cmp, no, nm))
  }
}

# filter on dN/dS alignment quality
dat_q <- dat[!is.na(dat$dN) & !is.na(dat$dS) & !is.na(dat$omega) &
               dat$mean_guidance_score >= min.score &
               dat$frac_cols_retained  >= min.frac, ]

# read each chromosome's GFF once
gff_cache <- setNames(
  lapply(c("juncea_X", "juncea_Y", "anax_12", "anax_X"), readGenes),
  c("juncea_X", "juncea_Y", "anax_12", "anax_X"))

# fusion breakpoint: the cut on the X that best separates Anax 12 (neo) from Anax X
# (ancestral) anchors; middle of the plateau of equally good cuts
estimate_breakpoint <- function(top = "juncea_X",
                                neo_cmp = "juncea_X__vs__anax_12",
                                anc_cmp = "juncea_X__vs__anax_X") {
  g <- gff_cache[[top]]
  pos <- function(cmp) {
    # unfiltered anchors (ratio_min = 0)
    sp <- synteny_partner(top, cmp, ratio_min = 0)
    if (is.null(sp) || !nrow(sp)) return(numeric(0))
    p <- g$mid[match(sp$top.gene, g$gene_id)]; p[!is.na(p)]
  }
  neo <- pos(neo_cmp); anc <- pos(anc_cmp)
  if (!length(neo) || !length(anc)) return(NULL)
  cuts <- sort(unique(c(neo, anc)))
  err  <- vapply(cuts, function(cx) mean(neo > cx) + mean(anc < cx), numeric(1))
  best <- cuts[err == min(err)]
  mean(range(best))
}
bp_est <- estimate_breakpoint()
if (!is.null(bp_est)) {
  ENT$juncea_X$bp <- bp_est
  message(sprintf("[breakpoint] using synteny-anchor estimate: %.3f Mb", bp_est / 1e6))
} else {
  message(sprintf("[breakpoint] estimator failed, falling back to hard-coded %.3f Mb",
                  ENT$juncea_X$bp / 1e6))
}

# breakpoint.txt is read by figs/fig4/fig4.R, figs/fig4/supp_gene_traffic.R and
# divergence_dating/02.select_loci.R. pairs.sbatch writes it too, so an existing file is
# verified (error if > bp_tol apart) and then used, never overwritten.
bp_file  <- file.path(run_dir, "breakpoint.txt")
bp_tol   <- 1000        # bp
if (!is.null(bp_est)) {
  bp_new <- as.integer(round(ENT$juncea_X$bp))
  if (file.exists(bp_file)) {
    bp_old <- as.integer(gsub("[^0-9]", "", readLines(bp_file)[1]))
    if (abs(bp_new - bp_old) > bp_tol)
      stop(sprintf("breakpoint disagreement: estimated %d, %s holds %d (%d bp apart). ",
                   bp_new, bp_file, bp_old, abs(bp_new - bp_old)),
           "One is stale. Delete the file to accept the estimate.")
    message(sprintf("[breakpoint] agrees with %s (%d bp, differ by %d); file left as is",
                    bp_file, bp_old, abs(bp_new - bp_old)))
    ENT$juncea_X$bp <- bp_old        # the file is authoritative once it exists
  } else {
    writeLines(sprintf("%d", bp_new), bp_file)
    message(sprintf("[breakpoint] wrote %s (%d bp)", bp_file, bp_new))
  }
} else {
  # never write the hard-coded fallback to breakpoint.txt
  warning("breakpoint estimator failed -- ", bp_file, " left untouched. ",
          "Fix the anchor/GFF inputs rather than accepting the fallback.")
  if (file.exists(bp_file))
    message(sprintf("[breakpoint] existing %s holds %s bp",
                    bp_file, readLines(bp_file)[1]))
}

# junction audit: X-Y gametologs on the ancestral X arm, with their Anax anchor
diagnose_junction <- function(bp = ENT$juncea_X$bp) {
  g   <- gff_cache[["juncea_X"]]
  xy  <- synteny_partner("juncea_X", "juncea_X__vs__juncea_Y")
  a12 <- synteny_partner("juncea_X", "juncea_X__vs__anax_12")$top.gene
  aX  <- synteny_partner("juncea_X", "juncea_X__vs__anax_X")$top.gene
  if (is.null(xy) || !nrow(xy)) return(invisible(NULL))
  xy$xmid <- g$mid[match(xy$top.gene, g$gene_id)]
  out <- xy[!is.na(xy$xmid) & xy$xmid >= bp, ]
  if (!nrow(out)) {
    message("[junction] no X-Y gametologs above ", round(bp/1e6, 2), " Mb")
    return(invisible(out))
  }
  out$anax_anchor <- ifelse(out$top.gene %in% a12, "Anax12 (neo)",
                     ifelse(out$top.gene %in% aX,  "AnaxX (ancestral)", "no Anax anchor"))
  out[order(out$xmid), c("top.gene", "gene_b", "xmid", "anax_anchor")]
}
print(diagnose_junction())

pal <- adjustcolor(rev(magma(101)), 0.6)

# Figure 1: Aeshna X vs Anax 12 / Anax X, MCScanX reciprocal 1:1 anchors
# (Anax 12 ribbons land on the neo-X arm, Anax X ribbons on the ancestral arm)
svg("Aeshna_X_vs_Anax_12_Anax_X.svg", width = 28, height = 7)
{
  fig1_title <- "Aeshna X  vs  Anax 12 / Anax X"
  fig1_top   <- "juncea_X"
  fig1_bottoms <- list(list(cmp = "juncea_X__vs__anax_12", align = "left"),
                       list(cmp = "juncea_X__vs__anax_X",  align = "bp"))
  topL <- fig1_top; topE <- ENT[[topL]]
  X_MB <- topE$len / 1e6

  # ribbon pairing: MCScanX reciprocal 1:1
  cmps_v <- vapply(fig1_bottoms, function(b) b$cmp, character(1))
  bp_tbl <- synteny_recip_partner(topL, cmps_v)

  # resolve each bottom into mapped 1:1 pairs + x-offset
  res <- list()
  for (b in fig1_bottoms) {
    parts  <- strsplit(b$cmp, "__vs__")[[1]]
    botL   <- setdiff(parts, topL)
    topcol <- if (topL == parts[1]) "gene_a" else "gene_b"
    botcol <- if (topL == parts[1]) "gene_b" else "gene_a"
    botE   <- ENT[[botL]]

    q <- if (!is.null(bp_tbl)) bp_tbl[bp_tbl$comparison == b$cmp, , drop = FALSE]
         else one2one(b$cmp)
    if (is.null(q$top.mid))
      q$top.mid <- gff_cache[[topL]]$mid[match(strip_tx(q[[topcol]]), gff_cache[[topL]]$gene_id)]
    if (is.null(q$bot.mid))
      q$bot.mid <- gff_cache[[botL]]$mid[match(strip_tx(q[[botcol]]), gff_cache[[botL]]$gene_id)]
    q <- q[!is.na(q$top.mid) & !is.na(q$bot.mid), ]

    al <- b$align
    off_mb <- if (is.null(al) || identical(al, "left")) 0 else
      if (identical(al, "right")) X_MB - botE$len / 1e6 else
        if (identical(al, "bp"))    topE$bp / 1e6 else al / 1e6
    res[[length(res) + 1]] <- list(E = botE, q = q, off_mb = off_mb)
  }

  xmax <- X_MB
  for (r in res) xmax <- max(xmax, r$off_mb + r$E$len / 1e6)

  par(mar = fig.mar)
  plot(NULL, xlim = c(0, xmax), ylim = c(-0.28, 1.25),
       xlab = "Position (Mb)", ylab = "", yaxt = "n", bty = "n", main = fig1_title)

  # bars
  draw_bar(topE, 1.0, 0)
  for (r in res) draw_bar(r$E, 0.0, r$off_mb)

  # left label
  text(-lab.frac * xmax, 1.0, topE$name, cex = 1.4, xpd = NA, adj = 1)

  # per-arm counts for each bottom
  bp_mb <- topE$bp / 1e6
  b1 <- res[[1]]; b2 <- res[[2]]
  dx <- 0.03 * xmax
  segs <- list(list(name = "neo-X",       cx = bp_mb / 2,          col = topE$col_neo, neo = TRUE),
               list(name = "ancestral X", cx = (bp_mb + X_MB) / 2, col = topE$col_anc, neo = FALSE))
  for (s in segs) {
    n1 <- sum(if (s$neo) b1$q$top.mid <  topE$bp else b1$q$top.mid >= topE$bp)
    n2 <- sum(if (s$neo) b2$q$top.mid <  topE$bp else b2$q$top.mid >= topE$bp)
    text(s$cx,      -0.13, s$name, col = s$col,    cex = 1.0, font = 2, xpd = NA)
    text(s$cx - dx, -0.21, n1,     col = b1$E$col, cex = 1.1, xpd = NA)
    text(s$cx + dx, -0.21, n2,     col = b2$E$col, cex = 1.1, xpd = NA)
  }

  # ribbons coloured by bottom chromosome, highest dS drawn first
  for (r in res) {
    d <- r$q[order(r$q$dS, decreasing = TRUE, na.last = FALSE), ]
    cols <- rep(adjustcolor(r$E$col, 0.5), nrow(d))
    for (i in seq_len(nrow(d))) {
      x0 <- d$top.mid[i] / 1e6
      x1 <- d$bot.mid[i] / 1e6 + r$off_mb
      bz <- bezierCubic(c(x0, 0.94), c(x0, 0.5), c(x1, 0.5), c(x1, 0.06))
      lines(bz$x, bz$y, col = cols[i], lwd = 0.5)
    }
  }
}
dev.off()

# Figure 2: Aeshna X vs Aeshna Y gametologs, MCScanX reciprocal 1:1 anchors,
# ribbons coloured by X arm (neo-X purple, ancestral X green)
svg("Aeshna_X_vs_Aeshna_Y.svg", width = 28, height = 7)
{
  fig2_title <- "Aeshna X  vs  Aeshna Y"
  fig2_top   <- "juncea_X"
  fig2_bottoms <- list(list(cmp = "juncea_X__vs__juncea_Y", align = "left"))
  topL <- fig2_top; topE <- ENT[[topL]]
  X_MB <- topE$len / 1e6

  # ribbon pairing: MCScanX reciprocal 1:1
  cmps_v <- vapply(fig2_bottoms, function(b) b$cmp, character(1))
  bp_tbl <- synteny_recip_partner(topL, cmps_v)

  # resolve each bottom into mapped 1:1 pairs + x-offset
  res <- list()
  for (b in fig2_bottoms) {
    parts  <- strsplit(b$cmp, "__vs__")[[1]]
    botL   <- setdiff(parts, topL)
    topcol <- if (topL == parts[1]) "gene_a" else "gene_b"
    botcol <- if (topL == parts[1]) "gene_b" else "gene_a"
    botE   <- ENT[[botL]]

    q <- if (!is.null(bp_tbl)) bp_tbl[bp_tbl$comparison == b$cmp, , drop = FALSE]
         else one2one(b$cmp)
    if (is.null(q$top.mid))
      q$top.mid <- gff_cache[[topL]]$mid[match(strip_tx(q[[topcol]]), gff_cache[[topL]]$gene_id)]
    if (is.null(q$bot.mid))
      q$bot.mid <- gff_cache[[botL]]$mid[match(strip_tx(q[[botcol]]), gff_cache[[botL]]$gene_id)]
    q <- q[!is.na(q$top.mid) & !is.na(q$bot.mid), ]

    al <- b$align
    off_mb <- if (is.null(al) || identical(al, "left")) 0 else
      if (identical(al, "right")) X_MB - botE$len / 1e6 else
        if (identical(al, "bp"))    topE$bp / 1e6 else al / 1e6
    res[[length(res) + 1]] <- list(E = botE, q = q, off_mb = off_mb)
  }

  xmax <- X_MB
  for (r in res) xmax <- max(xmax, r$off_mb + r$E$len / 1e6)

  par(mar = fig.mar)
  plot(NULL, xlim = c(0, xmax), ylim = c(-0.28, 1.25),
       xlab = "Position (Mb)", ylab = "", yaxt = "n", bty = "n", main = fig2_title)

  # bars and left label
  draw_bar(topE, 1.0, 0)
  for (r in res) draw_bar(r$E, 0.0, r$off_mb)

  text(-lab.frac * xmax, 1.0, topE$name, cex = 1.4, xpd = NA, adj = 1)

  # labels: neo-X vs ancestral X counts
  for (r in res) {
    cx <- r$off_mb + (r$E$len / 1e6) / 2
    n_neo <- sum(r$q$top.mid <  topE$bp)
    n_anc <- sum(r$q$top.mid >= topE$bp)
    qtr   <- (r$E$len / 1e6) / 4
    text(cx - qtr, -0.18, sprintf("neo-X\n(n = %d)", n_neo),
         col = topE$col_neo, cex = 1.0, xpd = NA)
    text(cx + qtr, -0.18, sprintf("ancestral X\n(n = %d)", n_anc),
         col = topE$col_anc, cex = 1.0, xpd = NA)
  }

  # ribbons coloured by top arm, highest dS drawn first
  for (r in res) {
    d <- r$q[order(r$q$dS, decreasing = TRUE, na.last = FALSE), ]
    cols <- adjustcolor(ifelse(d$top.mid < topE$bp, topE$col_neo, topE$col_anc), 0.5)
    for (i in seq_len(nrow(d))) {
      x0 <- d$top.mid[i] / 1e6
      x1 <- d$bot.mid[i] / 1e6 + r$off_mb
      bz <- bezierCubic(c(x0, 0.94), c(x0, 0.5), c(x1, 0.5), c(x1, 0.06))
      lines(bz$x, bz$y, col = cols[i], lwd = 2)
    }
  }
}
dev.off()
