#!/usr/bin/env Rscript
# Summarises kmer.sh output per window class: auto, neoX_arm / ancX_arm (the X split at NEOX_BREAK;
# NA = one 'neoX' class), neoY, micro. Tables go to OUTDIR/tables, plots to OUTDIR/tables/report.pdf.
# usage: summarize_kmer_sex.R OUTDIR W MICRO NEOX MIN_CHR_LEN NPERM [NEOX_BREAK] [NEOX_ARM] [NEOY]
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 6) stop("usage: summarize_kmer_sex.R OUTDIR W MICRO NEOX MIN_CHR_LEN NPERM [NEOX_BREAK] [NEOX_ARM]")
out <- a[1]; W <- as.numeric(a[2]); MICRO <- a[3]; NEOX <- a[4]; MINLEN <- as.numeric(a[5]); B <- as.integer(a[6])
BREAK <- if (length(a) >= 7) suppressWarnings(as.numeric(a[7])) else NA_real_
ARM   <- if (length(a) >= 8 && nzchar(a[8])) a[8] else "left"
NEOY  <- if (length(a) >= 9) a[9] else ""
set.seed(1)
tb <- function(f) file.path(out, "tables", f)
al <- function(f) file.path(out, "aln", f)
MIN_DEN <- 0.1 * W   # drop windows with fewer read-supported k-mer starts than this (gaps, repeat deserts)
CLS_ORDER <- c("auto", "ancX_arm", "neoX_arm", "neoX", "neoY", "micro")   # row order of the class tables

# chromosome classes on hap1: from step 0 metadata (chrom_classes.tsv), else by length
g <- fread(tb("hap1.genome"), header = FALSE, col.names = c("chrom", "len"))
if (nzchar(MICRO) && !MICRO %in% g$chrom) stop("MICRO not found in hap1: ", MICRO)
if (nzchar(NEOX) && !NEOX %in% g$chrom) warning("NEOX not found in hap1: ", NEOX)
if (nzchar(NEOY) && !NEOY %in% g$chrom) warning("NEOY not found in hap1: ", NEOY)
if (file.exists(tb("chrom_classes.tsv"))) {
  cls <- fread(tb("chrom_classes.tsv"), colClasses = list(character = c("seqid", "chr_name", "seq_name")))
  g <- merge(g, cls[, .(chrom = seqid, chr_name, class)], by = "chrom", all.x = TRUE)
  g[is.na(class), class := "unplaced"]
} else {
  g[, `:=`(chr_name = NA_character_,
           class = fifelse(chrom == MICRO, "micro",
                           fifelse(chrom == NEOX, "neoX", fifelse(len >= MINLEN, "auto", "unplaced"))))]
  cat("no chrom_classes.tsv - classifying by MIN_CHR_LEN =", MINLEN, "\n")
}
# the MICRO / NEOX / NEOY arguments always win over chrom_classes.tsv
if (nzchar(MICRO)) g[chrom == MICRO, class := "micro"]
if (nzchar(NEOX))  g[chrom == NEOX,  class := "neoX"]
if (nzchar(NEOY))  g[chrom == NEOY,  class := "neoY"]
g[, label := fifelse(is.na(chr_name) | chr_name == "", chrom, chr_name)]
if (g[class == "excluded", .N]) cat("excluded from all statistics:", g[class == "excluded", paste(label, collapse = ", ")], "\n")

# split the X into arms only when a breakpoint is given
SPLIT <- nzchar(NEOX) && NEOX %in% g$chrom && !is.na(BREAK)
if (SPLIT) cat(sprintf("X split at %.3f Mb on %s; neo-X arm is the %s side\n", BREAK / 1e6, NEOX, ARM))
if (nzchar(NEOX) && !SPLIT) cat("X treated as a single block (no breakpoint given)\n")

# per-window class: the chromosome's class, except on the X, which is split by position
add_wclass <- function(dt) {
  dt <- merge(dt, g[, .(chrom, class)], by = "chrom")
  dt[, wclass := class]
  if (SPLIT) dt[chrom == NEOX, wclass := fifelse((start < BREAK) == (ARM == "left"), "neoX_arm", "ancX_arm")]
  dt[]
}
# groups for per-chromosome tables and plots (the X contributes two groups when split), labels and colours
grp_levels <- c(g[class == "auto"][order(-len), chrom],
                if (SPLIT) c("neoX_arm", "ancX_arm") else if (nzchar(NEOX)) intersect(NEOX, g$chrom),
                if (nzchar(NEOY)) intersect(NEOY, g$chrom),
                if (nzchar(MICRO)) MICRO)
grp_label <- function(v) {
  m <- g[match(v, chrom), label]
  fifelse(is.na(m), v, m)
}
grp_class <- function(v) {
  cl <- g[match(v, chrom), class]
  fifelse(is.na(cl), v, cl)
}
grp_labels <- grp_label(grp_levels)
grp_cols   <- c(auto = "grey40", neoX = "dodgerblue3", neoX_arm = "dodgerblue3",
                ancX_arm = "darkorchid3", neoY = "darkorange2", micro = "firebrick",
                unplaced = "grey80", excluded = "grey90")[grp_class(grp_levels)]
add_grp <- function(dt) dt[, grp := fifelse(chrom == NEOX & SPLIT, wclass, chrom)][]

# chromosome order and offsets along the genome-wide plots (unique: no duplicate offset joins)
chr_plot <- unique(c(g[class == "auto"][order(-len), chrom],
                     if (nzchar(NEOX)) NEOX, if (nzchar(NEOY)) NEOY, if (nzchar(MICRO)) MICRO))
gk <- g[match(chr_plot, chrom)]
gk[, offset := cumsum(as.numeric(shift(len, fill = 0)))]   # as.numeric: integer cumsum overflows past 2.1 Gb

# male-specific proportion per window (num / den), single-copy track and, if present, all-k-mer track
read_win <- function(f, val) fread(tb(f), header = FALSE, col.names = c("chrom", "start", "end", val))
load_track <- function(num, den, label) {
  x <- merge(read_win(den, "den"), read_win(num, "num"), by = c("chrom", "start", "end"), all.x = TRUE)
  x[is.na(num), num := 0]
  x[, `:=`(prop = num / den, track = label)]
  x
}
tr <- load_track("windows_num1.tsv", "windows_den1.tsv", "single_copy")
if (file.exists(tb("windows_numA.tsv"))) tr <- rbind(tr, load_track("windows_numA.tsv", "windows_denA.tsv", "all_kmers"))
tr <- add_grp(add_wclass(tr))[den >= MIN_DEN & !wclass %in% c("unplaced", "excluded")]

# per-group and per-class summaries, fold vs the autosomal median
auto_med <- tr[wclass == "auto", .(auto_median = median(prop)), by = track]
chr_sum <- tr[, .(n_win = .N, median_prop = median(prop), mean_prop = mean(prop)), by = .(track, grp, wclass)]
chr_sum <- merge(chr_sum, auto_med, by = "track")
chr_sum[, `:=`(label = grp_label(grp), fold_vs_auto = median_prop / auto_median)]
chr_sum <- chr_sum[order(track, match(grp, grp_levels))]

class_sum <- tr[, .(n_win = .N, median_prop = median(prop), pooled_prop = sum(num) / sum(den)), by = .(track, wclass)]
class_sum <- merge(class_sum, auto_med, by = "track")
class_sum[, fold_vs_auto := median_prop / auto_median]
class_sum <- class_sum[order(track, match(wclass, CLS_ORDER))]

# block permutation per track and class: how often a random contiguous autosomal block of the same
# number of windows has a median >= the class's observed median (one-sided)
targets <- intersect(c("micro", "neoY", "neoX_arm", "ancX_arm", "neoX"), unique(tr$wclass))
perm <- list()
for (tk in unique(tr$track)) {
  x <- tr[track == tk]
  for (cl in targets) {
    obs <- x[wclass == cl][order(chrom, start), prop]; n <- length(obs)
    obs_med <- NA_real_; perm_med <- NA_real_; p_blk <- NA_real_; note <- "too few windows"
    if (n >= 3) {
      obs_med <- median(obs)
      note <- "no autosome long enough for a block of this size"
      auto <- x[wclass == "auto"][order(chrom, start)]
      pos <- auto[, .(nstart = max(0L, .N - n + 1L)), by = chrom][nstart > 0]
      if (nrow(pos) > 0) {
        # every possible block start (chromosome, offset), sampled B times
        props <- split(auto$prop, auto$chrom)
        starts <- rep(pos$chrom, pos$nstart)
        off <- sequence(pos$nstart)
        pick <- sample.int(length(starts), B, replace = TRUE)
        null_med <- numeric(length(pick))
        for (j in seq_along(pick)) {
          i <- pick[j]
          null_med[j] <- median(props[[starts[i]]][off[i]:(off[i] + n - 1)])
        }
        perm_med <- median(null_med)
        p_blk <- (1 + sum(null_med >= obs_med)) / (B + 1)
        note <- sprintf("%d possible autosomal blocks of %d windows", length(starts), n)
      }
    }
    perm[[length(perm) + 1]] <- data.table(track = tk, wclass = cl, n_win = n, obs_median = obs_med,
                                           perm_median = perm_med, p_block = p_blk, note = note)
  }
}
perm <- rbindlist(perm)

# composition control: proportion and depth ratio depend on how repetitive a window is, so everything
# is also reported within strata of single-copy content (single_copy / all_kmers denominators)
SCBRK <- c(-Inf, 0.2, 0.4, 0.6, Inf); SCLAB <- c("<0.2", "0.2-0.4", "0.4-0.6", ">0.6")
have_sc <- all(c("single_copy", "all_kmers") %in% tr$track)
if (have_sc) {
  scw <- dcast(tr[, .(chrom, start, track, den)], chrom + start ~ track, value.var = "den")
  scw <- scw[!is.na(single_copy) & !is.na(all_kmers) & all_kmers > 0]
  scw[, sc := single_copy / all_kmers][, scbin := cut(sc, SCBRK, labels = SCLAB)]
  scw <- scw[, .(chrom, start, sc, scbin)]
  trs <- merge(tr, scw, by = c("chrom", "start"))

  # median proportion per stratum, and fold vs the autosomes in the same stratum
  strat <- trs[, .(median_prop = median(prop), n = .N), by = .(track, wclass, scbin)]
  strat <- merge(strat, strat[wclass == "auto", .(track, scbin, auto_prop = median_prop)], by = c("track", "scbin"))
  strat[, fold_in_stratum := median_prop / auto_prop]
  strat <- strat[order(track, match(wclass, CLS_ORDER), scbin)]

  # direct standardisation: autosomal stratum medians weighted by the class's stratum composition
  std <- strat[, .(n_win = sum(n), observed = weighted.mean(median_prop, n),
                   expected_from_composition = weighted.mean(auto_prop, n)), by = .(track, wclass)]
  std[, fold_adjusted := observed / expected_from_composition]
  std <- merge(std, class_sum[, .(track, wclass, fold_unadjusted = fold_vs_auto)], by = c("track", "wclass"))
  std <- std[order(track, match(wclass, CLS_ORDER))]

  # composition-matched resampling: autosomal windows drawn to match the class's stratum composition
  pm <- list()
  for (tk in unique(trs$track)) {
    x <- trs[track == tk]
    for (cl in targets) {
      o <- x[wclass == cl]
      p_matched <- NA_real_
      if (nrow(o) >= 10) {
        tab <- o[, .N, by = scbin]
        pool <- split(x[wclass == "auto", prop], x[wclass == "auto", scbin])
        obs <- median(o$prop)
        Bm <- min(B, 10000)
        draws <- replicate(Bm, {
          vals <- numeric(0)
          for (i in seq_len(nrow(tab))) {
            pv <- pool[[as.character(tab$scbin[i])]]
            if (!is.null(pv) && length(pv)) vals <- c(vals, sample(pv, tab$N[i], replace = TRUE))
          }
          median(vals)
        })
        p_matched <- (1 + sum(draws >= obs)) / (Bm + 1)
      }
      pm[[length(pm) + 1]] <- data.table(track = tk, wclass = cl, p_matched = p_matched)
    }
  }
  pm <- rbindlist(pm)
  perm <- merge(perm, pm, by = c("track", "wclass"), all.x = TRUE)
  setnames(perm, "p_block", "p_block_unadjusted", skip_absent = TRUE)
}

# read depth per window, normalised to the autosomal median, per sex and by group / class
load_depth <- function(prefix) {
  f <- al(paste0(prefix, "_hap1.regions.bed.gz"))
  # fread() needs R.utils to open .gz directly; pipe through gzip instead so the dependency is not required
  d <- fread(cmd = paste("gzip -cd", shQuote(f)), header = FALSE, col.names = c("chrom", "start", "end", "depth"))
  d <- add_grp(add_wclass(d))[!wclass %in% c("unplaced", "excluded")]
  d[, norm := depth / median(d[wclass == "auto", depth])]
  d[]
}
have_depth <- file.exists(al("male_hap1.regions.bed.gz")) && file.exists(al("fem_hap1.regions.bed.gz"))
if (have_depth) {
  dm <- load_depth("male"); df <- load_depth("fem")
  depth_chr <- merge(dm[, .(male_norm = median(norm)), by = .(grp, wclass)],
                     df[, .(fem_norm = median(norm)), by = .(grp, wclass)], by = c("grp", "wclass"))
  depth_chr[, `:=`(label = grp_label(grp), fem_over_male = fem_norm / male_norm)]
  depth_chr <- depth_chr[order(match(grp, grp_levels))]
  depth_cls <- merge(dm[, .(male_norm = median(norm)), by = wclass],
                     df[, .(fem_norm = median(norm)), by = wclass], by = "wclass")
  depth_cls[, fem_over_male := fem_norm / male_norm]
  depth_cls <- depth_cls[order(match(wclass, CLS_ORDER))]
}

# female:male depth ratio within strata of single-copy content, and composition-adjusted
depth_strat <- NULL; depth_std <- NULL
if (have_depth && have_sc) {
  dw <- merge(dm[, .(chrom, start, wclass, male = norm)], df[, .(chrom, start, fem = norm)], by = c("chrom", "start"))
  dw <- merge(dw, scw, by = c("chrom", "start"))[wclass != "excluded"]
  dw[, ratio := fem / pmax(male, 1e-3)]
  depth_strat <- dw[, .(median_ratio = median(ratio), median_male = median(male), median_fem = median(fem), n = .N),
                    by = .(wclass, scbin)]
  depth_strat <- merge(depth_strat, depth_strat[wclass == "auto", .(scbin, auto_ratio = median_ratio)], by = "scbin")
  depth_strat <- depth_strat[order(match(wclass, CLS_ORDER), scbin)]
  depth_std <- depth_strat[, .(n_win = sum(n), observed = weighted.mean(median_ratio, n),
                               expected_from_composition = weighted.mean(auto_ratio, n)), by = wclass]
  depth_std[, fold_adjusted := observed / expected_from_composition]
}

# reciprocal control: female-specific k-mers on female contigs, each contig assigned to a hap1
# chromosome by primary alignment (>= 50 % of its aligned bases to one chromosome)
have_recip <- file.exists(al("fem_to_hap1.slim.tsv")) && file.exists(tb("windows_numF.tsv"))
PAFCOLS <- c("qname", "qlen", "qstart", "qend", "strand", "tname", "tlen", "tstart", "tend", "matches", "alnlen", "mapq", "tp")
if (have_recip) {
  paf <- fread(al("fem_to_hap1.slim.tsv"), header = FALSE, col.names = PAFCOLS)
  asg <- paf[tp == "P", .(bases = sum(matches)), by = .(qname, tname)][order(qname, -bases)]
  asg <- asg[, .(hap1_chrom = tname[1], frac_to_best = bases[1] / sum(bases), aligned_bases = sum(bases)), by = qname]
  asg[frac_to_best < 0.5, hap1_chrom := "ambiguous"]
  fw <- merge(read_win("windows_denF.tsv", "den"), read_win("windows_numF.tsv", "num"), by = c("chrom", "start", "end"), all.x = TRUE)
  fw[is.na(num), num := 0]
  fc <- fw[, .(num = sum(num), den = sum(den)), by = .(contig = chrom)]
  fc <- merge(fc, asg[, .(contig = qname, hap1_chrom)], by = "contig", all.x = TRUE)
  fc[is.na(hap1_chrom), hap1_chrom := "unaligned"]
  fc[, prop := num / den]
  fc <- merge(fc, g[, .(hap1_chrom = chrom, class)], by = "hap1_chrom", all.x = TRUE)
  recip <- fc[den >= MIN_DEN, .(n_contigs = .N, median_contig_prop = median(prop), pooled_prop = sum(num) / sum(den), kmers = sum(den)),
              by = .(hap1_chrom, class)]
  recip[, auto_pooled := recip[class == "auto", weighted.mean(pooled_prop, kmers)]]
  recip[, `:=`(label = grp_label(hap1_chrom), fold_vs_auto = pooled_prop / auto_pooled)]
  recip <- recip[order(match(hap1_chrom, grp_levels))]

  # contig-level Wilcoxon test of each non-autosomal assignment vs autosome-assigned contigs
  fcc <- fc[den >= MIN_DEN]
  auto_props <- fcc[class == "auto", prop]
  recip_test <- list()
  for (h in setdiff(unique(fcc$hap1_chrom), fcc[class == "auto", unique(hap1_chrom)])) {
    x <- fcc[hap1_chrom == h, prop]
    recip_test[[length(recip_test) + 1]] <- data.table(
      hap1_chrom = h, label = grp_label(h), n_contigs = length(x), median_contig_prop = median(x),
      auto_median = median(auto_props),
      p_wilcox_greater = if (length(x) >= 3) wilcox.test(x, auto_props, alternative = "greater")$p.value else NA_real_)
  }
  recip_test <- rbindlist(recip_test)

  # W candidates: female contigs with the highest female-specific proportion (size-filtered)
  w_cand <- fcc[order(-prop)][1:min(20, .N), .(contig, hap1_chrom, label = grp_label(hap1_chrom), class, kmers = den, prop)]
}

# where each candidate (micro, neoY) aligns in the rest of hap1, and which X arm it hits
# (a neo-Y should concentrate on the neo-X arm, above the fraction expected at random)
arm_tables <- list(); top_tables <- list()
for (role in c("micro", "neoy")) {
  f <- al(paste0(role, "_to_hap1rest.slim.tsv"))
  if (!file.exists(f)) next
  mp <- fread(f, header = FALSE, col.names = PAFCOLS)
  if (!nrow(mp)) next
  tp <- mp[, .(aligned_bases = sum(matches), n_aln = .N), by = tname][order(-aligned_bases)]
  tp[, `:=`(label = grp_label(tname), class = grp_class(tname), role = role)]
  top_tables[[role]] <- head(tp, 8)
  if (!nzchar(NEOX) || !SPLIT) next
  mx <- mp[tname == NEOX]
  if (!nrow(mx)) next
  mx[, arm := fifelse((tstart < BREAK) == (ARM == "left"), "neoX_arm", "ancX_arm")]
  at <- mx[, .(n_aln = .N, aligned_bases = sum(matches), median_identity = median(matches / alnlen)), by = arm]
  xlen <- g[chrom == NEOX, len]
  la <- if (ARM == "left") BREAK else xlen - BREAK
  at[, `:=`(role = role,
            expected_frac_if_random = fifelse(arm == "neoX_arm", la / xlen, 1 - la / xlen),
            observed_frac = aligned_bases / sum(aligned_bases))]
  arm_tables[[role]] <- at
}
micro_arm <- if (length(arm_tables)) rbindlist(arm_tables) else NULL
mp_top    <- if (length(top_tables)) rbindlist(top_tables) else NULL

# heterozygosity (step 9) per window and class, per 1000 callable sites in the same depth band
# (~0 in a hemizygous chromosome, autosomal level in a mosaic one)
have_het <- file.exists(tb("male_het_windows.tsv")) && file.exists(tb("fem_het_windows.tsv")) &&
  file.exists(tb("het_regions.txt")) && have_depth
if (have_het) {
  regs <- fread(tb("het_regions.txt"), header = FALSE)$V1
  # full window grid for the called regions (windows with zero hets are absent from the tables)
  grid <- list()
  for (ch in regs) {
    L <- g[chrom == ch, len]; n <- floor(L / W)
    if (n >= 1) grid[[length(grid) + 1]] <- data.table(chrom = ch, start = (0:(n - 1)) * W, end = (1:n) * W)
  }
  grid <- rbindlist(grid)

  # het counts and callable sites per window, with the window's normalised depth, per sex
  MIN_CALLABLE <- 1000   # windows with too few in-band sites give an unstable ratio
  het <- list()
  for (sx in c("male", "fem")) {
    dd <- if (sx == "male") dm else df
    h <- fread(tb(paste0(sx, "_het_windows.tsv")), header = FALSE,
               col.names = c("chrom", "start", "end", "n_het", "mean_altfrac", "n_balanced"))
    x <- merge(grid, h, by = c("chrom", "start", "end"), all.x = TRUE)
    x[is.na(n_het), `:=`(n_het = 0L, n_balanced = 0L)]
    cf <- tb(paste0(sx, "_callable_windows.tsv"))
    if (file.exists(cf)) {
      cal <- fread(cf, header = FALSE, col.names = c("chrom", "start", "end", "n_callable"))
      x <- merge(x, cal, by = c("chrom", "start", "end"), all.x = TRUE)
      x[is.na(n_callable), n_callable := 0L]
    } else x[, n_callable := NA_integer_]
    x <- merge(x, dd[, .(chrom, start, norm)], by = c("chrom", "start"))
    x <- add_grp(add_wclass(x))[!wclass %in% c("unplaced", "excluded")]
    x[, `:=`(sex = sx, het_per_mb = n_het / (end - start) * 1e6,
             het_per_kb_callable = fifelse(n_callable > 0, n_het / n_callable * 1e3, NA_real_))]
    het[[sx]] <- x
  }
  het <- rbindlist(het)

  # per class, relative to the autosomes of the same sex
  het_cls <- het[, .(n_win = .N, n_callable = sum(n_callable), n_het = sum(n_het),
                     het_per_mb = sum(n_het) / sum(end - start) * 1e6,
                     het_per_kb_callable = sum(n_het) / max(1, sum(n_callable)) * 1e3,
                     frac_balanced = sum(n_balanced) / max(1, sum(n_het)),
                     mean_altfrac = weighted.mean(mean_altfrac, n_het, na.rm = TRUE)), by = .(sex, wclass)]
  het_cls <- merge(het_cls, het_cls[wclass == "auto", .(sex, auto_het = het_per_kb_callable)], by = "sex")
  het_cls[, ratio_vs_auto := het_per_kb_callable / auto_het]
  het_cls <- het_cls[order(sex, match(wclass, CLS_ORDER))]
  # classes with too few callable sites in the depth band to trust the ratio
  het_thin <- het_cls[n_callable < 5e5, .(sex, wclass, n_callable)]
}

# assembly-to-assembly coverage of each hap1 chromosome
cov_files <- c(fem_primary = "asmcov_fem_on_hap1_primary.tsv", fem_all = "asmcov_fem_on_hap1_all.tsv",
               hap2_primary = "asmcov_hap2_on_hap1_primary.tsv")
asmcov <- g[, .(chrom, class, len)]
for (n in names(cov_files)) {
  if (file.exists(tb(cov_files[n]))) {
    x <- fread(tb(cov_files[n]), header = FALSE, col.names = c("chrom", "len2", "covered", "frac"))
    asmcov <- merge(asmcov, x[, .(chrom, frac)], by = "chrom", all.x = TRUE)
    setnames(asmcov, "frac", paste0("frac_cov_", n))
  }
}
asmcov <- asmcov[!class %in% c("unplaced", "excluded")][order(match(chrom, gk$chrom))]
asmcov[, label := grp_label(chrom)]

# coverage of each candidate by the rest of hap1 and by hap2
q_files <- c(micro_vs_rest_of_hap1 = "asmcov_micro_by_hap1rest_all.tsv", micro_vs_hap2 = "asmcov_micro_by_hap2_all.tsv",
             neoy_vs_rest_of_hap1  = "asmcov_neoy_by_hap1rest_all.tsv",  neoy_vs_hap2  = "asmcov_neoy_by_hap2_all.tsv")
micro_q <- list()
for (cmp in names(q_files)) {
  if (file.exists(tb(q_files[[cmp]]))) {
    micro_q[[cmp]] <- fread(tb(q_files[[cmp]]), header = FALSE, col.names = c("chrom", "len", "covered", "frac"))
  }
}
micro_q <- rbindlist(micro_q, idcol = "comparison")

# light contamination screen: GC per chromosome vs the autosomes
if (file.exists(tb("hap1_seqstats.tsv"))) {
  ss <- fread(tb("hap1_seqstats.tsv"), header = FALSE, col.names = c("chrom", "len", "gc"))
  ss <- merge(ss, g[, .(chrom, class, label)], by = "chrom")
  amu <- ss[class == "auto", mean(gc)]; asd <- ss[class == "auto", sd(gc)]
  ss[, gc_z_vs_auto := (gc - amu) / asd]
  ss[, flag := fifelse(abs(gc_z_vs_auto) > 3, "GC outlier", "")]
  seq_flags <- ss[!class %in% c("unplaced", "excluded")][order(match(chrom, gk$chrom))]
}

# write tables
fwrite(chr_sum, tb("chrom_summary.tsv"), sep = "\t")
fwrite(class_sum, tb("class_summary.tsv"), sep = "\t")
if (have_sc) {
  fwrite(strat, tb("kmer_by_stratum.tsv"), sep = "\t")
  fwrite(std, tb("kmer_composition_adjusted.tsv"), sep = "\t")
}
if (!is.null(depth_strat)) {
  fwrite(depth_strat, tb("depth_by_stratum.tsv"), sep = "\t")
  fwrite(depth_std, tb("depth_composition_adjusted.tsv"), sep = "\t")
}
fwrite(tr[order(track, match(grp, grp_levels), start)], tb("windows_male_specific.tsv"), sep = "\t")
fwrite(perm, tb("permutation.tsv"), sep = "\t")
if (have_depth) {
  fwrite(depth_chr, tb("depth_by_chrom.tsv"), sep = "\t")
  fwrite(depth_cls, tb("depth_by_class.tsv"), sep = "\t")
}
if (have_recip) {
  fwrite(recip, tb("female_specific_by_chrom.tsv"), sep = "\t")
  fwrite(fc, tb("fem_contig_assignment.tsv"), sep = "\t")
  fwrite(recip_test, tb("female_specific_test.tsv"), sep = "\t")
  fwrite(w_cand, tb("female_specific_top_contigs.tsv"), sep = "\t")
}
if (have_het) {
  fwrite(het_cls, tb("het_by_class.tsv"), sep = "\t")
  fwrite(het, tb("het_windows.tsv"), sep = "\t")
}
fwrite(asmcov, tb("asm_coverage_by_chrom.tsv"), sep = "\t")
if (!is.null(micro_arm)) fwrite(micro_arm, tb("micro_alignment_arms.tsv"), sep = "\t")
if (exists("seq_flags")) fwrite(seq_flags, tb("seq_flags.tsv"), sep = "\t")

# console summary
options(width = 200, datatable.print.nrows = 200)
if (file.exists(tb("set_sizes.tsv"))) {
  cat("\n== set sizes / dropout ==\n")
  # integer64="double": k-mer counts exceed .Machine$integer.max and print as garbage unless bit64 is loaded
  print(fread(tb("set_sizes.tsv"), integer64 = "double")); print(fread(tb("dropout.tsv"), integer64 = "double"))
}
cat("\n== male-specific proportion by class ==\n"); print(class_sum)
if (have_sc) {
  cat("\n== COMPOSITION-ADJUSTED male-specific enrichment ==\n")
  cat("   fold_unadjusted compares a class to the autosomal median regardless of how repetitive it is.\n")
  cat("   The autosomal proportion itself varies several-fold across single-copy strata, so a class that is\n")
  cat("   more repetitive than the autosomes is inflated by composition alone. fold_adjusted applies the\n")
  cat("   class's own stratum composition to the autosomal stratum medians: THIS is the number to quote.\n")
  print(std[, .(track, wclass, n_win, observed, expected_from_composition, fold_adjusted, fold_unadjusted)])
  cat("\n== male-specific proportion within strata of single-copy content ==\n")
  print(dcast(strat, track + wclass ~ scbin, value.var = "fold_in_stratum"))
  cat("   (window counts per stratum)\n")
  print(dcast(strat, track + wclass ~ scbin, value.var = "n"))
}

cat("\n== per-chromosome (X shown as two arms) ==\n")
print(chr_sum[wclass != "auto" | grp %in% grp_levels[1:min(5, length(grp_levels))],
              .(track, label, wclass, n_win, median_prop, auto_median, fold_vs_auto)])
cat("\n== permutation vs autosomal windows ==\n")
cat("   p_block_unadjusted: contiguous autosomal blocks - controls window autocorrelation, NOT composition.\n")
if (have_sc) cat("   p_matched: autosomal windows resampled to the class's single-copy composition - the reverse.\n")
print(perm)
if (have_depth) {
  cat("\n== read depth / autosomal median (male / female) ==\n")
  cat("   auto 1.0/1.0 | ancX_arm 0.5/1.0 | neoX_arm 0.5-1.0/1.0 | male-limited 0.5/~0 |\n")
  cat("   X-like (hemizygous in male) 0.5/1.0 | Z-like (hemizygous in female) 1.0/0.5 | mosaic: intermediate, uniform along the chromosome\n")
  print(depth_cls); print(depth_chr[, .(label, wclass, male_norm, fem_norm, fem_over_male)])
}
if (!is.null(depth_strat)) {
  cat("\n== COMPOSITION-ADJUSTED female:male depth ratio ==\n")
  cat("   The autosomal ratio itself falls sharply in repeat-rich windows (a mapping/MAPQ artifact present on\n")
  cat("   every chromosome), so a repeat-rich class looks female-depleted for that reason alone.\n")
  print(depth_std)
  cat("\n== female:male depth ratio within strata of single-copy content ==\n")
  print(dcast(depth_strat, wclass ~ scbin, value.var = "median_ratio"))
  cat("   (window counts per stratum)\n")
  print(dcast(depth_strat, wclass ~ scbin, value.var = "n"))
}
if (have_recip) {
  cat("\n== female-specific direction: proportion by assigned hap1 chromosome ==\n")
  cat("   'unaligned' and 'ambiguous' are bins to LOOK AT, not discard: a diverged W-like homolog lands there\n")
  print(recip[, .(label, class, n_contigs, pooled_prop, auto_pooled, fold_vs_auto)])
  cat("\n== female-specific: contig-level Wilcoxon vs autosome-assigned contigs ==\n")
  print(recip_test)
  cat("\n== female contigs with the highest female-specific proportion (W candidates if they cluster) ==\n")
  print(w_cand)
}
if (have_het) {
  cat("\n== heterozygosity by class ==\n")
  cat("   het_per_kb_callable is the comparable statistic: het SNPs per 1000 sites inside the depth band,\n")
  cat("   so power is matched across classes of different depth. ratio_vs_auto ~0 = hemizygous, ~1 = diploid/mosaic.\n")
  print(het_cls[, .(sex, wclass, n_callable, n_het, het_per_kb_callable, ratio_vs_auto, het_per_mb, frac_balanced, mean_altfrac)])
  if (nrow(het_thin)) {
    cat("\n   CAUTION - too few callable sites in the depth band to trust these rows:\n")
    print(het_thin)
  }
}
if (!is.null(mp_top)) {
  cat("\n== where each candidate aligns in the rest of hap1 ==\n")
  print(mp_top[, .(role, label, class, n_aln, aligned_bases)])
}
if (!is.null(micro_arm)) {
  cat("\n== candidate alignments to the X, by arm (a true neo-Y should hit the neo-X arm) ==\n")
  print(micro_arm[, .(role, arm, n_aln, aligned_bases, median_identity, observed_frac, expected_frac_if_random)])
}
cat("\n== assembly coverage of hap1 chromosomes ==\n"); print(asmcov[, -c("class", "len"), with = FALSE])
if (nrow(micro_q)) {
  cat("\n== micro covered by diverged homology elsewhere ==\n")
  print(micro_q)
}
if (exists("seq_flags")) {
  cat("\n== GC per chromosome ==\n")
  print(seq_flags[, .(label, len, gc, gc_z_vs_auto, flag)])
}

# plots: genome-wide male-specific proportion
pdf(tb("report.pdf"), width = 11, height = 7)
xbreak <- if (SPLIT) gk[chrom == NEOX, offset] + BREAK else NA_real_
t1 <- merge(tr[track == "single_copy"], gk[, .(chrom, offset)], by = "chrom")
t1[, x := offset + start]
plot(t1$x, t1$prop, pch = 16, cex = 0.4, col = grp_cols[match(t1$grp, grp_levels)], xaxt = "n",
     xlab = "hap1 (autosomes by size, then X, then micro)", ylab = "fraction of single-copy male k-mers unseen in female",
     main = sprintf("Male-specific k-mer proportion, %d kb windows", W / 1000))
abline(v = gk$offset, col = "grey85")
if (SPLIT) abline(v = xbreak, col = "darkorchid3", lty = 3, lwd = 2)
abline(h = auto_med[track == "single_copy", auto_median], lty = 2)
axis(1, at = gk$offset + gk$len / 2, labels = gk$label, las = 2, cex.axis = 0.6)

# per-group window distributions, one page per track
for (tk in unique(tr$track)) {
  boxplot(prop ~ factor(grp, levels = grp_levels, labels = grp_labels), data = tr[track == tk], las = 2, cex.axis = 0.6,
          col = grp_cols, xlab = "", ylab = "male-specific proportion",
          main = sprintf("Per-chromosome window distributions (%s); X split into arms", tk))
  abline(h = auto_med[track == tk, auto_median], lty = 2)
}

# genome-wide depth of both sexes, and normalised depth per group
if (have_depth) {
  d2 <- merge(rbind(dm[, .(chrom, start, norm, sex = "male")], df[, .(chrom, start, norm, sex = "female")]),
              gk[, .(chrom, offset)], by = "chrom")
  d2[, x := offset + start]
  plot(d2[sex == "male", x], d2[sex == "male", norm], pch = 16, cex = 0.3, col = "grey30", ylim = c(0, 2), xaxt = "n",
       xlab = "hap1", ylab = "window depth / autosomal median", main = "Read depth on hap1 (grey = male, red = female)")
  points(d2[sex == "female", x], d2[sex == "female", norm], pch = 16, cex = 0.3, col = "firebrick")
  abline(v = gk$offset, col = "grey85"); abline(h = c(0.5, 1), lty = 2)
  if (SPLIT) abline(v = xbreak, col = "darkorchid3", lty = 3, lwd = 2)
  axis(1, at = gk$offset + gk$len / 2, labels = gk$label, las = 2, cex.axis = 0.6)
  m <- t(as.matrix(depth_chr[match(grp_levels, grp), .(male_norm, fem_norm)])); colnames(m) <- grp_labels
  barplot(m, beside = TRUE, las = 2, cex.names = 0.6, col = c("grey30", "firebrick"), legend.text = c("male", "female"),
          ylab = "median window depth / autosomal median", main = "Normalised depth (X split into arms)")
  abline(h = c(0.5, 1), lty = 2)
}

# reciprocal control per assigned hap1 chromosome
if (have_recip) {
  r2 <- recip[match(intersect(grp_levels, recip$hap1_chrom), hap1_chrom)]
  barplot(r2$pooled_prop, names.arg = r2$label, las = 2, cex.names = 0.6,
          col = grp_cols[match(r2$hap1_chrom, grp_levels)],
          ylab = "female-specific proportion (pooled over assigned contigs)", main = "Reciprocal control")
  abline(h = recip[, unique(auto_pooled)], lty = 2)
}

# composition-controlled enrichment per track, and depth ratio by stratum
if (have_sc) {
  par(mfrow = c(1, 2))
  for (tk in unique(strat$track)) {
    st <- strat[track == tk & wclass != "auto"]
    if (!nrow(st)) next
    m <- dcast(st, wclass ~ scbin, value.var = "fold_in_stratum")
    mm <- as.matrix(m[, -1, with = FALSE]); rownames(mm) <- m$wclass
    barplot(mm, beside = TRUE, las = 1, col = grp_cols[match(rownames(mm), grp_class(grp_levels))],
            legend.text = rownames(mm), args.legend = list(x = "topleft", bty = "n", cex = 0.7),
            xlab = "single-copy content of window", ylab = "fold vs autosomes IN THE SAME STRATUM",
            main = sprintf("Composition-controlled enrichment (%s)", tk))
    abline(h = 1, lty = 2)
  }
  if (!is.null(depth_strat)) {
    ds <- depth_strat[wclass != "auto"]
    m <- dcast(ds, wclass ~ scbin, value.var = "median_ratio")
    mm <- as.matrix(m[, -1, with = FALSE]); rownames(mm) <- m$wclass
    barplot(mm, beside = TRUE, las = 1, col = grp_cols[match(rownames(mm), grp_class(grp_levels))],
            legend.text = rownames(mm), args.legend = list(x = "topleft", bty = "n", cex = 0.7),
            xlab = "single-copy content of window", ylab = "median female:male depth ratio",
            main = "Depth ratio by composition")
    aline <- depth_strat[wclass == "auto"][order(scbin), median_ratio]
    points(x = (seq_along(aline) - 1) * (nrow(mm) + 1) + (nrow(mm) + 1) / 2, y = aline, pch = 4, cex = 1.4, lwd = 2)
    legend("topright", pch = 4, legend = "autosomes", bty = "n", cex = 0.7)
  }
  par(mfrow = c(1, 1))
}

# heterozygosity per group, one panel per sex
if (have_het) {
  par(mfrow = c(1, 2))
  for (sx in c("male", "fem")) {
    hs <- het[sex == sx]; lv <- intersect(grp_levels, unique(hs$grp))
    boxplot(het_per_kb_callable ~ factor(grp, levels = lv, labels = grp_label(lv)), data = hs[n_callable >= MIN_CALLABLE],
            las = 2, cex.axis = 0.6, col = grp_cols[match(lv, grp_levels)], xlab = "",
            ylab = "het SNPs per 1000 callable sites", main = sprintf("Heterozygosity, %s reads on hap1", sx))
    abline(h = het_cls[sex == sx & wclass == "auto", het_per_kb_callable], lty = 2)
  }
  par(mfrow = c(1, 1))
}

# assembly-to-assembly coverage per chromosome
if (nrow(asmcov) && any(grepl("frac_cov_", names(asmcov)))) {
  m <- t(as.matrix(asmcov[, grep("frac_cov_", names(asmcov)), with = FALSE])); colnames(m) <- asmcov$label
  barplot(m, beside = TRUE, las = 2, cex.names = 0.6, col = c("firebrick", "salmon", "steelblue")[seq_len(nrow(m))],
          legend.text = sub("frac_cov_", "", rownames(m)), ylab = "fraction of hap1 chromosome covered",
          main = "Assembly-to-assembly coverage")
}
invisible(dev.off())
cat("\nwrote", tb("report.pdf"), "\n")
