# Fig. S2: depth, CDS depth and repeat panels along the microchromosome, one svg each
library(data.table)

ROOT     <- "../.."
RUNKMER  <- "microchromosome/sex_kmers/kmer_sex_out"  # aln/{male,fem}_hap1.regions.bed.gz, tables/windows_den1.tsv
RUNDEPTH <- "microchromosome/cds_depth"               # mchrom_k1_depth_per_gene.tsv
CDS_BED  <- "microchromosome/cds_depth/cds.annevo.bed"
RUNREP   <- "microchromosome/repeats"
REP_TSV  <- "repeats_chroms.tsv.gz"
OUTDIR   <- "."
STEMS    <- c(depth = "figS2a_depth", cds = "figS2b_cds", repeats = "figS2c_repeats")

MICRO <- "OZ369553.1"; TNAME <- "microchromosome"
AUTO  <- c("OZ369542.1","OZ369543.1","OZ369544.1","OZ369545.1","OZ369546.1","OZ369547.1",
           "OZ369548.1","OZ369549.1","OZ369550.1","OZ369551.1","OZ369552.1")

# running median and band (depth panel)
SMOOTH_WIN  <- 20                # 100 kb windows (20 = 2 Mb)
BAND_Q      <- c(0.025, 0.975)
EXP_K       <- 101               # autosomal windows nearest in k-mer count behind each expected value
EPS         <- 0.02              # pseudocount on normalised depth before log2
MIN_RUN_WIN <- 5                 # shortest run below the band that is listed (windows)

# CDS
CDS_K <- 1; MIN_BP <- 300; CAP_CDS <- 3; PT_CEX_CDS <- 1.1

# repeats
REP_WIN <- 100000
KIM_BRK <- c(0, 5, 10, 15, 20, 30, 100)
KIM_LAB <- c("0-5", "5-10", "10-15", "15-20", "20-30", ">30")
REP_TOP_LINES <- 3.4

# appearance
PANEL_W  <- 4700; PANEL_H <- 1600
SCALE    <- PANEL_H / 760         # text, ticks, points, lines, margins relative to 1900 x 760 panels (or a number)
TICK_LEN <- -0.5                  # tick length in text lines; more negative = longer
TOP_LINES <- 2.2                 # legend margin for the depth and CDS panels
SEX_COL  <- c(fem = "#8da0cb", male = "#fc8d62")
PT_CEX   <- 0.55; PT_ALPHA <- 0.6
PT_CEX_BG <- 0.45; PT_ALPHA_BG <- 0.35   # background points under a running median
CAP_DEPTH <- 2
set.seed(1)

# paths under ROOT, input check, plain or gzipped tables
p <- function(...) { parts <- c(path.expand(ROOT), ...); do.call(file.path, as.list(parts[nzchar(parts)])) }
need <- function(f) { if (!file.exists(f)) stop("not found: ", f, "\n  working directory: ", getwd(),
                                                "\n  check ROOT and the paths at the top of this script"); f }
read_tab <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("gzip -cd", shQuote(f)), ...) else fread(f, ...)

# svg device; pointsize scales text, margins, ticks and symbols (line widths are scaled where drawn)
dev_any <- function(f, w, h) svg(f, width = w / 190, height = h / 190, pointsize = 12 * SCALE)
LINE_PX <- local({                                   # px per margin line on the device in use
  tf <- tempfile(); dev_any(tf, PANEL_W, PANEL_H)
  par(mar = c(0, 0, 1, 0)); v <- par("mai")[3]; invisible(dev.off()); unlink(tf); v * 190
})

# open a panel: same plotting area in every panel, the top margin adds to the height
open_panel <- function(stem, top) {
  f <- file.path(OUTDIR, paste0(stem, ".svg"))
  dev_any(f, PANEL_W, PANEL_H + round((top - 0.8) * LINE_PX))
  par(mar = c(4.0, 4.8, top, 0.8), mgp = c(2.9, 0.7, 0), tcl = TICK_LEN, lwd = SCALE); f   # lwd: box, default lines
}
ax <- function(...) axis(..., lwd = SCALE, lwd.ticks = SCALE)          # axis line and tick marks, scaled
darken <- function(col, f = 0.6) { v <- col2rgb(col) * f; rgb(v[1], v[2], v[3], maxColorValue = 255) }
rmed <- function(v, k) frollapply(v, k, median, align = "center")
top_legend <- function(...) { u <- par("usr")
  legend(x = mean(u[1:2]), y = u[4], xjust = 0.5, yjust = 0, xpd = NA, horiz = TRUE, bty = "n", cex = 0.8, ...) }

# window depth per sex, normalised to the autosomal median
dep <- list()
for (sx in c("male", "fem")) {
  d <- read_tab(need(p(RUNKMER, "aln", paste0(sx, "_hap1.regions.bed.gz"))),
                col.names = c("chrom", "start", "end", "depth"))[chrom %in% c(AUTO, MICRO)]
  dep[[sx]] <- d[, norm := depth / median(d[chrom %in% AUTO, depth])][]
}
D <- merge(dep$male[, .(chrom, start, end, male = norm)], dep$fem[, .(chrom, start, fem = norm)],
           by = c("chrom", "start"))

# add single-copy k-mer counts (windows without any are dropped), log scales
den <- fread(need(p(RUNKMER, "tables", "windows_den1.tsv")), col.names = c("chrom", "start", "end", "den"))
D <- merge(D, den[, .(chrom, start, den)], by = c("chrom", "start"))[den > 0]
setorder(D, chrom, start)
D[, `:=`(lf = log2(fem + EPS), lm = log2(male + EPS), lx = log10(den + 1))]
T12  <- D[chrom == MICRO]
LEN  <- max(T12$end); LENM <- LEN / 1e6

# expected log depth from k-mer count (running median over the EXP_K autosomal windows nearest in count)
for (v in c("lf", "lm")) {
  a <- D[chrom %in% AUTO & !is.na(get(v))][order(lx)]
  cur <- a[, sm := runmed(get(v), EXP_K, endrule = "median")][, .(y = median(sm)), by = lx]
  D[, (paste0("e_", v)) := approx(cur$lx, cur$y, xout = lx, rule = 2)$y]
  D[, (paste0("s_", v)) := get(v) - get(paste0("e_", v))]
}

# per sex: running medians of depth and of its expectation; the band is the BAND_Q quantiles of
# autosomal (curve - expected curve) added to the expected curve
S <- list()
for (sx in c("fem", "male")) {
  v <- c(fem = "lf", male = "lm")[[sx]]
  X <- D[!is.na(get(v)), .(chrom, start, end, val = get(v), e = get(paste0("e_", v)), s = get(paste0("s_", v)))]
  X[, `:=`(c = rmed(val, SMOOTH_WIN), ce = rmed(e, SMOOTH_WIN)), by = chrom]
  X[, cs := c - ce]
  q <- X[chrom %in% AUTO, quantile(cs, BAND_Q, na.rm = TRUE)]
  X[, `:=`(lo = ce + q[1], hi = ce + q[2])][, below := c < lo]
  S[[sx]] <- list(target = X[chrom == MICRO], auto_below = X[chrom %in% AUTO, mean(below, na.rm = TRUE)])
}

# runs of >= MIN_RUN_WIN windows below the band; the female's are written out
runs_below <- function(m) {
  r <- rle(m$below %in% TRUE); e <- cumsum(r$lengths); s <- e - r$lengths + 1
  data.table(s, e, v = r$values)[v == TRUE & (e - s + 1) >= MIN_RUN_WIN,
                                 .(from_mb = m$start[s] / 1e6, to_mb = m$end[e] / 1e6, windows = e - s + 1)]
}
REG <- runs_below(S$fem$target)
fwrite(REG, file.path(OUTDIR, paste0(STEMS[["depth"]], "_regions.tsv")), sep = "\t")

# x axis in Mb and window midpoints, shared by the depth and CDS panels
xaxis <- function() {
  at <- pretty(c(0, LENM), n = 10); at <- at[at <= LENM]
  ax(1, at = at, cex.axis = 0.9); mtext(sprintf("position on %s (Mb)", TNAME), side = 1, line = 2.5, cex = 0.95)
}
xm <- function(d) (d$start + d$end) / 2e6
written <- character()

# depth panel: window points (capped) behind the running medians of both sexes
f <- open_panel(STEMS[["depth"]], TOP_LINES)
plot(NA, xlim = c(0, LENM), ylim = c(0, CAP_DEPTH * 1.04), axes = FALSE, xaxs = "i", yaxs = "i",
     xlab = "", ylab = sprintf("normalized depth, %.1f Mb running median", SMOOTH_WIN / 10))
pts <- rbind(T12[, .(x = xm(T12), y = fem, s = "fem")], T12[, .(x = xm(T12), y = male, s = "male")])[sample(.N)]
points(pts$x, pmin(pts$y, CAP_DEPTH), pch = 16, cex = PT_CEX_BG, col = adjustcolor(SEX_COL[pts$s], PT_ALPHA_BG))
for (sx in c("fem", "male")) {
  m <- S[[sx]]$target; ok <- !is.na(m$c)
  lines(xm(m)[ok], pmin(CAP_DEPTH, 2^m$c[ok]), col = darken(SEX_COL[[sx]]), lwd = 2.4 * SCALE)
}
ax(2, at = c(0, .5, 1, 1.5, 2), las = 1, cex.axis = 0.9); xaxis(); box(bty = "l")
L <- data.table(lab = c("female", "male"), lty = 1, lwd = 2.4 * SCALE, pch = NA_real_, cex = NA_real_,
                col = c(darken(SEX_COL[["fem"]]), darken(SEX_COL[["male"]])))
top_legend(legend = L$lab, lty = L$lty, lwd = L$lwd, pch = L$pch, pt.cex = L$cex, col = L$col,
           text.width = strwidth(L$lab, cex = 0.8) * 1.15, x.intersp = 0.8)
invisible(dev.off()); written <- c(written, f)

# CDS panel: per-gene depth at the shared-ortholog genes, normalised to the autosomal median
if (file.exists(p(CDS_BED))) {
  G <- fread(need(p(RUNDEPTH, sprintf("mchrom_k%d_depth_per_gene.tsv", CDS_K))))[cds_bp >= MIN_BP &
                                                                                  chrom %in% c(AUTO, MICRO)]
  G <- G[shared == 1]
  bed <- fread(p(CDS_BED), header = FALSE, select = 1:4, col.names = c("chrom", "start", "end", "gene"))
  bed[, gene := sub("\\.t\\d+$", "", gene)]
  G <- merge(G, bed[, .(start = min(start), end = max(end)), by = gene], by = "gene")
  G[, `:=`(male_n = male / median(G[chrom %in% AUTO, male]), fem_n = fem / median(G[chrom %in% AUTO, fem]))]
  g <- G[chrom == MICRO]
  f <- open_panel(STEMS[["cds"]], TOP_LINES)
  plot(NA, xlim = c(0, LENM), ylim = c(0, CAP_CDS * 1.04), axes = FALSE, xaxs = "i", yaxs = "i",
       xlab = "", ylab = "normalized CDS depth")
  pts <- rbind(g[, .(x = xm(g), y = fem_n, s = "fem")], g[, .(x = xm(g), y = male_n, s = "male")])[sample(.N)]
  points(pts$x, pmin(pts$y, CAP_CDS), pch = 16, cex = PT_CEX_CDS, col = adjustcolor(SEX_COL[pts$s], PT_ALPHA))
  ax(2, at = 0:3, las = 1, cex.axis = 0.9); xaxis(); box(bty = "l")
  top_legend(legend = c("female", "male"), pch = 16, pt.cex = 1.3, col = SEX_COL)
  invisible(dev.off()); written <- c(written, f)
  cat(sprintf("[cds] %d genes on the %s, normalised to %d autosomal genes of the same set\n",
              nrow(g), TNAME, G[chrom %in% AUTO, .N]))
} else message("note: CDS_BED not found - the CDS panel is skipped")

# repeat panel: proportion of each window masked, by Kimura divergence bin
fr <- p(RUNREP, REP_TSV)
if (file.exists(fr)) {
  REP <- read_tab(fr, na.strings = "NA")[chrom == MICRO]
  REP[, kbin := factor(as.character(cut(kimura, KIM_BRK, labels = KIM_LAB, right = FALSE)), levels = c(KIM_LAB, "NA"))]
  REP[is.na(kbin), kbin := "NA"]

  # each base goes to ONE hit (highest score if present, else longest), so a window cannot exceed 1
  pr <- if ("score" %in% names(REP)) as.numeric(REP$score) else as.numeric(REP$end - REP$start + 1)
  H  <- data.table(s0 = as.numeric(REP$start - 1), e1 = as.numeric(REP$end) - 1, pr = pr, kbin = REP$kbin)
  br <- sort(unique(c(H$s0, H$e1 + 1, seq(0, LEN, by = REP_WIN), LEN)))
  SEG <- data.table(s = br[-length(br)], e = br[-1] - 1)
  setkey(H, s0, e1)
  ov  <- foverlaps(SEG, H, by.x = c("s", "e"), by.y = c("s0", "e1"), type = "within", nomatch = 0L)
  win <- ov[ov[, .I[which.max(pr)], by = .(s, e)]$V1][, .(bp = e - s + 1, w = s %/% REP_WIN, kbin)]

  # bp per bin and window, as a fraction of window length
  nwin <- ceiling(LEN / REP_WIN); wlen <- pmin(REP_WIN, LEN - (0:(nwin - 1)) * REP_WIN)
  levs <- c(KIM_LAB, "NA")
  Mk <- matrix(0, length(levs), nwin, dimnames = list(levs, NULL))
  pc <- win[, .(bp = sum(bp)), by = .(kbin, w)]
  for (i in seq_len(nrow(pc)))
    Mk[as.character(pc$kbin[i]), pc$w[i] + 1] <- Mk[as.character(pc$kbin[i]), pc$w[i] + 1] + pc$bp[i]
  Mk <- sweep(Mk, 2, wlen, "/")
  KIM_COL <- local({
    v <- if (requireNamespace("viridisLite", quietly = TRUE)) viridisLite::viridis(length(KIM_LAB))
         else hcl.colors(length(KIM_LAB), "viridis")
    c(rev(v), "#BBBBBB")                             # young repeats bright
  })

  # stacked bars, axes in Mb, legend of the bins present
  f <- open_panel(STEMS[["repeats"]], REP_TOP_LINES)
  barplot(Mk, space = 0, border = NA, col = KIM_COL, axes = FALSE, names.arg = rep("", ncol(Mk)),
          xlab = "", ylab = sprintf("repeat proportion in %gkb windows", REP_WIN / 1e3),
          ylim = c(0, 1), xaxs = "i", yaxs = "i")
  ax(2, at = seq(0, 1, 0.25), las = 1, cex.axis = 0.9)
  at <- pretty(c(0, LENM), n = 10); at <- at[at <= LENM]
  ax(1, at = at * 1e6 / REP_WIN, labels = at, cex.axis = 0.9)
  mtext(sprintf("position on %s (Mb)", TNAME), side = 1, line = 2.5, cex = 0.95)
  box(bty = "l")
  keep <- rowSums(Mk) > 0
  top_legend(legend = levs[keep], fill = KIM_COL[keep], border = NA, title = "Kimura divergence (%)",
             x.intersp = 0.6, text.width = max(strwidth(levs[keep], cex = 0.8)))
  invisible(dev.off()); written <- c(written, f)
} else message("note: ", fr, " not found - the repeat panel is skipped")

# summary: share of the microchromosome below the band (autosomal share = false-positive rate), runs below
cat(sprintf("\nrunning median %.1f Mb; band = %g-%g%% of repeat-matched autosomal sequence\n",
            SMOOTH_WIN / 10, 100 * BAND_Q[1], 100 * BAND_Q[2]))
cat(sprintf("  female below her band: %.1f%% of the %s   (autosomes: %.1f%% -- the false-positive rate)\n",
            100 * S$fem$target[, mean(below, na.rm = TRUE)], TNAME, 100 * S$fem$auto_below))
cat(sprintf("  male below his band:   %.1f%%   (autosomes: %.1f%%)\n",
            100 * S$male$target[, mean(below, na.rm = TRUE)], 100 * S$male$auto_below))
cat(sprintf("\ngrey regions (female below band, >= %d windows):\n", MIN_RUN_WIN)); print(REG)
RM <- runs_below(S$male$target)
cat(sprintf("\nmale below his band (not shaded), >= %d windows:\n", MIN_RUN_WIN)); print(RM)
fwrite(RM, file.path(OUTDIR, paste0(STEMS[["depth"]], "_regions_male.tsv")), sep = "\t")
cat("\nwrote", written, file.path(OUTDIR, paste0(STEMS[["depth"]], "_regions.tsv")), sep = "\n  "); cat("\n")
