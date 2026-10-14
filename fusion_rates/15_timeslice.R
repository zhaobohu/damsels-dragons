#!/usr/bin/env Rscript
# Time-slice test: is the focal clade enriched for sex-autosome fusions relative to every other
# lineage alive at its crown age T (and at multiples of T)? Parallel over maps (fork, L'Ecuyer-CMRG).
#   Rscript 15_timeslice.R [n_maps] [n_null] [cores]    (defaults: all, 9999, THREADS)

suppressPackageStartupMessages({ library(ape); library(data.table); library(parallel) })
source("config.R"); source("tree_utils.R")

args   <- commandArgs(trailingOnly = TRUE)
N_USE  <- if (length(args) >= 1) as.integer(args[1]) else NA_integer_
N_NULL <- if (length(args) >= 2) as.integer(args[2]) else 9999L
CORES  <- if (length(args) >= 3) as.integer(args[3]) else as.integer(THREADS)
# slice times as multiples of the focal crown age (only >= 1, so the focal clade stays whole)
T_FACTORS <- c(1.0, 1.05, 1.1, 1.25, 1.5, 2.0)
MAPS   <- "simmaps.rds"
OUT_TSV <- "timeslice_clades.tsv"; OUT_SUM <- "timeslice_summary.tsv"; OUT_FIG <- "fig_timeslice.pdf"
OUT_MAP <- "timeslice_permap.tsv"      # one row per map per T: rank, p_marg, p_maxT

tx <- fread(TAXA_TSV); tx[, tip := gsub(" ", "_", species)]
focal_tips <- gsub(" ", "_", FOCAL_TIPS)

# one map: slice clades and the relabelling null at every T factor
process_map <- function(i) {
  setDTthreads(1L)                        # forked workers must not each spawn threads
  m <- maps[[i]]; m$tip.label <- sub("^_R_", "", m$tip.label)
  if (!is.rooted(m)) return(NULL)

  # fusion events with their position along the branch (time since the parent node)
  ev_rows <- list()
  for (e in seq_along(m$maps)) {
    seg <- m$maps[[e]]; g <- names(seg); if (length(g) < 2) next
    cum <- cumsum(seg)
    for (k in seq_len(length(g) - 1L)) {
      fr <- strsplit(g[k], "_", fixed = TRUE)[[1]]; to <- strsplit(g[k + 1L], "_", fixed = TRUE)[[1]]
      if (length(fr) < 2 || length(to) < 2) next
      nf <- suppressWarnings(as.numeric(fr[1])); nt <- suppressWarnings(as.numeric(to[1]))
      if (is.na(nf) || is.na(nt) || nt >= nf) next          # fusions only
      ev_rows[[length(ev_rows) + 1L]] <- list(edge = e, pos = unname(cum[k]),
                                              is_SA = as.integer(fr[2] != to[2]), nchr = nf,
                                              sys = if (fr[2] == "X0") "XO" else "XY")
    }
  }
  if (!length(ev_rows)) return(NULL)
  ev <- rbindlist(ev_rows)
  ev[, p := p_null(nchr, sys)]

  # node ages, clade membership, focal crown age
  dep <- node.depth.edgelength(m); age <- max(dep) - dep
  dt <- tips_by_node(m); anc_list <- ancestors_by_node(m)
  ft <- intersect(focal_tips, m$tip.label); if (length(ft) < 2) return(NULL)
  T_focal <- age[getMRCA(m, ft)]

  # one null draw for every event (N_NULL x events), reused across T factors
  pvec <- ev$p
  lab  <- matrix(rbinom(N_NULL * length(pvec), 1L, rep(pvec, each = N_NULL)), nrow = N_NULL)

  out_cl <- list(); out_pm <- list()
  par_ <- m$edge[, 1]; chd <- m$edge[, 2]
  for (tf in T_FACTORS) {
    Tslice <- T_focal * tf

    # one slice clade per edge alive at Tslice (parent older, child at or younger)
    crossing <- which(age[par_] > Tslice & age[chd] <= Tslice)
    if (!length(crossing)) next
    evT <- copy(ev)
    evT[, age_ev := age[par_[edge]] - pos]                 # age at which the event happened
    evT[, below := age_ev <= Tslice]

    # per slice clade: events below Tslice on its stem plus all events inside it
    cl_rows <- list()
    for (e in crossing) {
      c0 <- chd[e]
      tips <- dt[[c0]]
      is_in <- logical(length(chd))
      for (k in seq_along(chd)) is_in[k] <- chd[k] != c0 && c0 %in% anc_list[[chd[k]]]
      inside <- which(is_in)
      ee <- evT[(edge == e & below) | edge %in% inside]
      S <- sum(ee$is_SA); K <- nrow(ee)
      E <- sum(ee$p); V <- sum(ee$p * (1 - ee$p))
      lmyr <- (Tslice - age[c0]) + sum(m$edge.length[inside])
      cl_rows[[length(cl_rows) + 1L]] <- data.table(
        edge = e, n_tips = length(tips),
        n_focal = sum(m$tip.label[tips] %in% focal_tips),
        S = S, A = K - S, K = K, E = E, V = V,
        z = if (V > 0) (S - E) / sqrt(V) else NA_real_,
        frac = if (K > 0) S / K else NA_real_,
        lmyr = lmyr, rate = if (lmyr > 0) S / lmyr else NA_real_,
        ev_idx = list(which((evT$edge == e & evT$below) | evT$edge %in% inside)))
    }
    cl <- rbindlist(cl_rows)

    # focal = the slice clade with the most focal tips; keep clades with K_MIN+ fusions and V > 0
    cl[, focal := n_focal == max(n_focal) & n_focal > 0]
    if (sum(cl$focal) > 1) cl[focal == TRUE, focal := seq_len(.N) == 1]
    keep <- cl[(K >= K_MIN | focal) & V > 0]
    if (!any(keep$focal) || nrow(keep) < 2) next
    fz <- keep[focal == TRUE, z]; if (is.na(fz)) next

    # null z for every clade: membership matrix (events x clades), S_null = lab %*% M
    M <- matrix(0, nrow = length(pvec), ncol = nrow(keep))
    for (j in seq_len(nrow(keep))) M[keep$ev_idx[[j]], j] <- 1
    Snull <- lab %*% M
    znull <- sweep(sweep(Snull, 2, keep$E, "-"), 2, sqrt(keep$V), "/")
    fcol  <- which(keep$focal)
    zmax  <- apply(znull, 1, max, na.rm = TRUE)
    hits_marg <- sum(znull[, fcol] >= fz); hits_max <- sum(zmax >= fz)

    # clade rows and the per-map row (ranks, p_marg, p_maxT)
    others <- keep[focal == FALSE]
    keep[, `:=`(iter = i, T_factor = tf, T_Ma = Tslice, ev_idx = NULL)]
    out_cl[[length(out_cl) + 1L]] <- keep
    out_pm[[length(out_pm) + 1L]] <- data.table(
      iter = i, T_factor = tf, T_Ma = Tslice, n_clades = nrow(keep),
      focal_tips = keep[focal == TRUE, n_tips], focal_in = keep[focal == TRUE, n_focal],
      focal_S = keep[focal == TRUE, S], focal_K = keep[focal == TRUE, K],
      focal_frac = keep[focal == TRUE, frac], focal_z = fz, focal_rate = keep[focal == TRUE, rate],
      rank_z    = rank(-keep$z, ties.method = "min")[keep$focal],
      rank_frac = rank(-keep$frac, ties.method = "min")[keep$focal],
      n_reach_z = sum(others$z >= fz, na.rm = TRUE),
      p_marg = (hits_marg + 1) / (N_NULL + 1), p_maxT = (hits_max + 1) / (N_NULL + 1))
  }
  list(cl = rbindlist(out_cl), pm = rbindlist(out_pm))
}

# run every map in parallel
maps <- readRDS(MAPS)
if (!is.na(N_USE)) maps <- maps[seq_len(min(N_USE, length(maps)))]
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
cat(sprintf("%d maps, %d null relabellings, T factors %s, K_MIN = %d, %d cores\n",
            length(maps), N_NULL, paste(T_FACTORS, collapse = "/"), K_MIN, CORES))
t0 <- Sys.time()
res <- mclapply(seq_along(maps), process_map, mc.cores = CORES, mc.set.seed = TRUE, mc.preschedule = TRUE)

# count worker errors, drop errored and empty maps, stack the tables
bad <- 0L; all_clades <- list(); per_map <- list()
for (r in res) {
  if (inherits(r, "try-error")) { bad <- bad + 1L; next }
  if (is.null(r)) next
  all_clades[[length(all_clades) + 1L]] <- r[["cl"]]
  per_map[[length(per_map) + 1L]] <- r[["pm"]]
}
if (bad > 0) cat(sprintf("warning: %d maps errored in workers\n", bad))
cat(sprintf("done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

C <- rbindlist(all_clades); P <- rbindlist(per_map)
fwrite(C, OUT_TSV, sep = "\t"); fwrite(P, OUT_MAP, sep = "\t")
if (!nrow(P)) stop("no usable maps", call. = FALSE)
missing_T <- setdiff(T_FACTORS, unique(P$T_factor))
if (length(missing_T))
  cat(sprintf("\nnote: T_factor %s produced no usable maps -- the slice fell above the root\n",
              paste(missing_T, collapse = ", ")))

# summary per T factor, medians across maps
q <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
S <- P[, .(maps = .N,
           T_Ma = round(median(T_Ma), 1),
           n_clades = as.numeric(median(n_clades)),
           focal_tips = as.numeric(median(focal_tips)), focal_in = as.numeric(median(focal_in)),
           focal_S = as.numeric(median(focal_S)), focal_K = as.numeric(median(focal_K)),
           frac = round(median(focal_frac), 3), frac_lo = round(q(focal_frac, .025), 3), frac_hi = round(q(focal_frac, .975), 3),
           z = round(median(focal_z), 2),
           rank_z = as.numeric(median(rank_z)), top_z = round(mean(rank_z == 1), 3),
           rank_frac = as.numeric(median(rank_frac)), top_frac = round(mean(rank_frac == 1), 3),
           clades_reach = round(mean(n_reach_z / pmax(n_clades - 1, 1)), 3),
           p_marg = round(median(p_marg), 4), p_maxT = round(median(p_maxT), 4),
           p_maxT_lo = round(q(p_maxT, .025), 4), p_maxT_hi = round(q(p_maxT, .975), 4)),
       by = T_factor][order(T_factor)]
fwrite(S, OUT_SUM, sep = "\t")

cat("\n=== time-slice comparison, median across maps ===\n")
cat("  T_factor 1.0 is the focal crown age; others are multiples of it.\n\n")
print(S[, .(T_factor, T_Ma, n_clades, focal_tips, focal_in, focal_S, focal_K, frac, z,
            rank_z, top_z, rank_frac, top_frac, p_marg, p_maxT)], row.names = FALSE)
cat("\n  n_clades     slice clades with >= K_MIN fusions (disjoint, so this IS the family size)\n")
cat("  focal_tips   tips in the focal slice clade; focal_in of them are focal species\n")
cat("  rank_z/frac  focal clade's rank by z / by SA fraction; top_* = fraction of maps ranked 1st\n")
cat("  p_maxT       max over disjoint clades -- 'corrected for n_clades' is literally true here\n")

# figure at T_factor 1: SA fraction vs fusions per slice clade (all maps, no jitter: stacked
# points saturate) over an approximate 95% band for random pairing
C1 <- C[T_factor == 1.0]
PT_CEX   <- 1.6      # point size
PT_ALPHA <- "12"     # hex fill alpha per point (~7%); ~40 stacked = saturated
if (nrow(C1)) {
  pdf(OUT_FIG, width = 7, height = 5)
  par(mar = c(4.5, 4.5, 2.5, 1))
  pbar <- C1[, sum(E) / sum(K)]                 # pooled null SA probability
  # normal approximation to the binomial interval (a guide band only; the test uses relabelling)
  ks <- seq(K_MIN, max(C1$K))
  se <- sqrt(pbar * (1 - pbar) / ks)
  lo <- pmax(0, pbar - 1.96 * se); hi <- pmin(1, pbar + 1.96 * se)
  plot(NA, xlim = range(C1$K), ylim = c(0, 1), log = "x",
       xlab = "fusions in the clade since T (all maps, T = focal crown age)",
       ylab = "sex-autosome fusions / all fusions",
       main = sprintf("%d maps; one point per slice clade per map", uniqueN(C1$iter)),
       cex.main = 0.9, font.main = 1)
  polygon(c(ks, rev(ks)), c(lo, rev(hi)), col = "#B4B2A955", border = NA)
  abline(h = pbar, col = "grey40", lty = 3)
  o <- C1[focal == FALSE]; f <- C1[focal == TRUE]
  points(o$K, o$frac, pch = 21, bg = paste0("#5F5E5A", PT_ALPHA), col = "#5F5E5A",
         cex = PT_CEX, lwd = 0.5)
  points(f$K, f$frac, pch = 21, bg = paste0("#D85A30", PT_ALPHA), col = "#993C1D",
         cex = PT_CEX * 1.15, lwd = 0.8)
  legend("topright",
         c("focal clade (one per map)", "other slice clades", "approx. 95% band, random pairing",
           "darker = more maps at that point"),
         pch = c(21, 21, 15, NA), pt.bg = c("#D85A30AA", "#5F5E5AAA", NA, NA),
         col = c("#993C1D", "#5F5E5A", "#B4B2A9", NA), pt.cex = c(1.5, 1.3, 1.8, NA),
         bty = "n", cex = 0.8)
  invisible(dev.off())
  cat(sprintf("\nwritten: %s, %s, %s, %s\n", OUT_TSV, OUT_MAP, OUT_SUM, OUT_FIG))
}
