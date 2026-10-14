# Anchor-partner rank filter for MCScanX anchors, shared by figs/fig3/fig3c.R,
# divergence_dating/02.select_loci.R and figs/fig4/fig4.R (change the rule here only).

suppressMessages(library(data.table))

# Rank statistics for every hit pair of one comparison, from comp/<cmp>/ab.rank.blast:
#   ratio    = bits(pair) / bits(query's best partner-chromosome hit)
#   rank_pos = position of the partner among the query's partner-chromosome hits (1 = best)
# Returns one row per ordered pair (both directions), or NULL if there is no search file.
rank_stats <- function(cmp, run_dir) {
  # uncapped search; fall back to ab.blast (-k 5, ranks truncated) if it is missing
  p <- file.path(run_dir, "comp", cmp, "ab.rank.blast")
  if (!file.exists(p)) {
    p2 <- file.path(run_dir, "comp", cmp, "ab.blast")
    if (!file.exists(p2)) {
      warning("[rank] no ab.rank.blast or ab.blast for ", cmp,
              " -- anchors cannot be rank-filtered")
      return(NULL)
    }
    message("[rank] ", cmp, ": no ab.rank.blast, falling back to ab.blast ",
            "(-k 5; ranks truncated at 5 hits per query)")
    p <- p2
  }
  parts <- strsplit(cmp, "__vs__")[[1]]
  ida <- file.path(run_dir, "chrom", paste0(parts[1], ".ids"))
  idb <- file.path(run_dir, "chrom", paste0(parts[2], ".ids"))

  # read hits (query, subject, pid, bitscore) and drop self hits
  bl <- fread(p, header = FALSE, select = c(1, 2, 3, 12))
  setnames(bl, c("gene_a", "gene_b", "pid", "bits"))
  bl <- bl[gene_a != gene_b]

  # keep cross-chromosome hits only, so same-chromosome paralogs stay out of the denominator
  if (file.exists(ida) && file.exists(idb)) {
    A <- readLines(ida); B <- readLines(idb)
    bl <- bl[(gene_a %chin% A & gene_b %chin% B) |
             (gene_a %chin% B & gene_b %chin% A)]
  } else {
    warning("[rank] missing .ids for ", cmp,
            " -- denominator includes same-chromosome paralogs")
  }
  if (!nrow(bl)) { warning("[rank] no cross-chromosome hits for ", cmp); return(NULL) }

  # best hit per pair, then ratio and rank within each true query
  bl <- bl[, .(pid = max(pid), bits = max(bits)), by = .(gene_a, gene_b)]
  bl[, ratio    := bits / max(bits), by = gene_a]
  bl[, rank_pos := frank(-bits, ties.method = "min"), by = gene_a]

  # merge both orientations on an unordered key (lenient: max ratio, min rank),
  # then emit both orderings so the caller's merge finds exactly one row
  bl[, `:=`(k1 = pmin(gene_a, gene_b), k2 = pmax(gene_a, gene_b))]
  agg <- bl[, .(ratio = max(ratio), pid = max(pid), rank_pos = min(rank_pos)),
            by = .(k1, k2)]
  rbind(agg[, .(gene_a = k1, gene_b = k2, ratio, pid, rank_pos)],
        agg[, .(gene_a = k2, gene_b = k1, ratio, pid, rank_pos)])
}

# Attach ratio/pid/rank_pos to anchors (transcript-level gene_a / gene_b) and drop failing ones.
# ratio_min = 0 attaches without filtering; pid_min > 0 also requires pid < pid_min to drop;
# na_keep = TRUE keeps pairs with no search row. Defaults come from min.ratio / min.pid /
# ratio.na.keep in the calling script if defined.
attach_rank <- function(d, cmp, run_dir,
                        ratio_min = if (exists("min.ratio")) min.ratio else 0.65,
                        pid_min   = if (exists("min.pid"))   min.pid   else 0,
                        na_keep   = if (exists("ratio.na.keep")) ratio.na.keep else FALSE) {
  if (is.null(d) || !nrow(d)) return(d)
  rs <- rank_stats(cmp, run_dir)
  if (is.null(rs)) {
    d$ratio <- NA_real_; d$pid <- NA_real_; d$rank_pos <- NA_integer_
    return(d)
  }
  was_dt <- is.data.table(d)
  d <- merge(as.data.frame(d), as.data.frame(rs),
             by = c("gene_a", "gene_b"), all.x = TRUE, sort = FALSE)

  # drop anchors with ratio below ratio_min (and pid below pid_min, if set)
  if (ratio_min > 0) {
    bad <- !is.na(d$ratio) & d$ratio < ratio_min
    if (pid_min > 0) bad <- bad & !is.na(d$pid) & d$pid < pid_min
    keep <- if (isTRUE(na_keep)) !bad else !bad & !is.na(d$ratio)
    d <- d[keep, , drop = FALSE]
  }
  if (was_dt) as.data.table(d) else d
}

# One-line summary of rank_pos for the methods text; returns the table of ranks (1, 2, 3+).
rank_report <- function(d, label = "anchors") {
  if (is.null(d) || !nrow(d) || is.null(d$rank_pos)) {
    cat("[rank report]", label, ": no rank_pos column\n"); return(invisible(NULL))
  }
  n <- nrow(d); r <- d$rank_pos
  cat(sprintf(
    "[rank report] %s: n = %d; partner is the query's best partner-chromosome hit for %d (%.1f%%), 2nd for %d, 3rd or worse for %d, no search row for %d\n",
    label, n, sum(r == 1, na.rm = TRUE), 100 * mean(r == 1, na.rm = TRUE),
    sum(r == 2, na.rm = TRUE), sum(r >= 3, na.rm = TRUE), sum(is.na(r))))
  invisible(table(pmin(r, 3), useNA = "ifany"))
}
