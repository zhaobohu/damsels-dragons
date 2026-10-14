# Fig. S4: repeat proportion in 100 kb windows by Kimura divergence, one svg per chromosome (X split at the junction)
library(data.table)

ROOT    <- "../.."
RUNREP  <- "microchromosome/repeats"                 # repeats_chroms.tsv.gz from run_repeats.sbatch
REP_TSV <- "repeats_chroms.tsv.gz"
RUNKMER <- "microchromosome/sex_kmers/kmer_sex_out"  # tables/hap1.genome (chromosome lengths)
OUTDIR  <- "."
WHICH   <- NULL                           # NULL = every region; or e.g. c("12", "X_neo", "Y")

AUTO  <- c("OZ369542.1","OZ369543.1","OZ369544.1","OZ369545.1","OZ369546.1","OZ369547.1",
           "OZ369548.1","OZ369549.1","OZ369550.1","OZ369551.1","OZ369552.1")
MICRO <- "OZ369553.1"; XCHR <- "OZ369554.1"; YCHR <- "OZ369555.1"
BP_FILE <- "gene_pairs_ds/breakpoint.txt"            # fusion junction on the X
X_COORDS <- "chromosome"                  # "chromosome": ancestral arm keeps its X coordinates
                                          # "arm": each arm's axis starts at 0

REP_WIN <- 100000
KIM_BRK <- c(0, 5, 10, 15, 20, 30, 100)
KIM_LAB <- c("0-5", "5-10", "10-15", "15-20", "20-30", ">30")
REP_TOP_LINES <- 3.4                      # top margin holding the legend; the panel grows by this margin

PANEL_W <- 2200; PANEL_H <- 860           # px at res 190, every panel the same size whatever its length

# paths under ROOT, input check, fusion junction, plain or gzipped tables
p <- function(...) { parts <- c(path.expand(ROOT), ...); do.call(file.path, as.list(parts[nzchar(parts)])) }
need <- function(f) { if (!file.exists(f)) stop("not found: ", f, "\n  working directory: ", getwd(),
                                                "\n  check ROOT and the paths at the top of this script"); f }
BP <- as.numeric(readLines(need(p(BP_FILE))))
read_tab <- function(f, ...) if (grepl("\\.gz$", f)) fread(cmd = paste("gzip -cd", shQuote(f)), ...) else fread(f, ...)
LINE_PX <- local({                        # px per margin line on the device in use
  tf <- tempfile()
  svg(tf, width = PANEL_W / 190, height = PANEL_H / 190)
  par(mar = c(0, 0, 1, 0)); v <- par("mai")[3]; invisible(dev.off()); unlink(tf); v * 190
})

# regions: autosomes 1-11, the microchromosome (12), the X arms either side of the junction, the neo-Y
glen <- fread(need(p(RUNKMER, "tables", "hap1.genome")), header = FALSE, col.names = c("chrom", "len"))
len_of <- function(ch) { v <- glen[chrom == ch, len]; if (!length(v)) stop(ch, " not in hap1.genome"); v }
REG <- rbind(
  data.table(name = as.character(seq_along(AUTO)), chrom = AUTO, from = 0, to = sapply(AUTO, len_of),
             xlab = sprintf("position on chromosome %d (Mb)", seq_along(AUTO))),
  data.table(name = "12", chrom = MICRO, from = 0, to = len_of(MICRO), xlab = "position on microchromosome (Mb)"),
  data.table(name = c("X_neo", "X_ancestral"), chrom = XCHR, from = c(0, BP), to = c(BP, len_of(XCHR)),
             xlab = c("position on X, neo-X arm (Mb)", "position on X, ancestral arm (Mb)")),
  data.table(name = "Y", chrom = YCHR, from = 0, to = len_of(YCHR), xlab = "position on neo-Y (Mb)"))
if (!is.null(WHICH)) {
  bad <- setdiff(WHICH, REG$name); if (length(bad)) stop("unknown region(s): ", paste(bad, collapse = ", "),
                                                        "; choose from ", paste(REG$name, collapse = ", "))
  REG <- REG[name %in% WHICH]
}

# repeat hits (needed columns and chromosomes only) with a Kimura bin and an overlap priority
fr   <- need(p(RUNREP, REP_TSV))
hdr  <- names(read_tab(fr, nrows = 0))
cols <- intersect(c("chrom", "start", "end", "kimura", "score"), hdr)
REP  <- read_tab(fr, select = cols, na.strings = "NA")[chrom %in% REG$chrom]
REP[, kbin := factor(as.character(cut(kimura, KIM_BRK, labels = KIM_LAB, right = FALSE)), levels = c(KIM_LAB, "NA"))]
REP[is.na(kbin), kbin := "NA"]
REP[, pr := if ("score" %in% cols) as.numeric(score) else as.numeric(end - start + 1)]   # priority uses the FULL hit
cat(sprintf("[repeats] %s hits on %d chromosomes; overlaps resolved by %s\n", format(nrow(REP), big.mark = ","),
            uniqueN(REP$chrom), if ("score" %in% cols) "alignment score" else "hit length"))

LEVS <- c(KIM_LAB, "NA")
KIM_COL <- local({
  v <- if (requireNamespace("viridisLite", quietly = TRUE)) viridisLite::viridis(length(KIM_LAB))
       else hcl.colors(length(KIM_LAB), "viridis")
  c(rev(v), "#BBBBBB")                     # young repeats bright
})

written <- character(); tabs <- list(); summ <- list()
for (i in seq_len(nrow(REG))) {
  r <- REG[i]
  if (!r$chrom %in% REP$chrom) { message("skipped ", r$name, ": no repeat hits for ", r$chrom, " in ", fr); next }

  # hits clipped to the region; windows of REP_WIN start at the region's start
  from <- r$from; to <- r$to
  R <- REP[chrom == r$chrom & end > from & start - 1 < to]
  nwin <- ceiling((to - from) / REP_WIN); wlen <- pmin(REP_WIN, (to - from) - (0:(nwin - 1)) * REP_WIN)
  M <- matrix(0, length(LEVS), nwin, dimnames = list(LEVS, NULL))
  if (nrow(R)) {
    H  <- R[, .(s0 = pmax(as.numeric(start) - 1, from), e1 = pmin(as.numeric(end), to) - 1, pr, kbin)][s0 <= e1]

    # elementary segments (every hit edge and window edge is a breakpoint), each given to its top-priority hit
    br <- sort(unique(c(H$s0, H$e1 + 1, seq(from, to, by = REP_WIN), to)))
    br <- br[br >= from & br <= to]
    SEG <- data.table(s = head(br, -1), e = br[-1] - 1)
    setkey(H, s0, e1)
    ov  <- foverlaps(SEG, H, by.x = c("s", "e"), by.y = c("s0", "e1"), type = "within", nomatch = 0L)
    win <- ov[ov[, .I[which.max(pr)], by = .(s, e)]$V1][, .(bp = e - s + 1, w = (s - from) %/% REP_WIN, kbin)]

    # bp per Kimura bin and window
    pc  <- win[, .(bp = sum(bp)), by = .(kbin, w)]
    for (j in seq_len(nrow(pc)))
      M[as.character(pc$kbin[j]), pc$w[j] + 1] <- M[as.character(pc$kbin[j]), pc$w[j] + 1] + pc$bp[j]
  }
  M <- sweep(M, 2, wlen, "/")
  if (any(colSums(M) > 1 + 1e-9)) warning(r$name, ": a window exceeds 1 after overlap resolution - should not happen")

  # draw: stacked bars, axes in Mb, legend of the bins present above the plot
  f <- file.path(OUTDIR, paste0(r$name, ".svg"))
  svg(f, width = PANEL_W / 190, height = (PANEL_H + round((REP_TOP_LINES - 0.8) * LINE_PX)) / 190)
  par(mar = c(4.0, 4.8, REP_TOP_LINES, 0.8), mgp = c(2.9, 0.7, 0))
  barplot(M, space = 0, border = NA, col = KIM_COL, axes = FALSE, names.arg = rep("", ncol(M)),
          xlab = "", ylab = sprintf("repeat proportion in %gkb windows", REP_WIN / 1e3),
          ylim = c(0, 1), xaxs = "i", yaxs = "i")
  axis(2, at = seq(0, 1, 0.25), las = 1, cex.axis = 0.9)
  off <- if (X_COORDS == "chromosome") r$from else 0            # value shown at the left edge
  lo <- off / 1e6; hi <- (off + r$to - r$from) / 1e6
  at <- pretty(c(lo, hi), n = 10); at <- at[at >= lo & at <= hi]
  axis(1, at = (at - lo) * 1e6 / REP_WIN, labels = at, cex.axis = 0.9)
  mtext(r$xlab, side = 1, line = 2.5, cex = 0.95)
  box(bty = "l")
  keep <- rowSums(M) > 0; u <- par("usr")
  legend(x = mean(u[1:2]), y = u[4], xjust = 0.5, yjust = 0, xpd = NA, horiz = TRUE, bty = "n", cex = 0.8,
         legend = LEVS[keep], fill = KIM_COL[keep], border = NA, title = "Kimura divergence (%)",
         x.intersp = 0.6, text.width = max(strwidth(LEVS[keep], cex = 0.8)))
  invisible(dev.off())
  written <- c(written, f)

  # per-window table, and the region's masked share overall and in the youngest/oldest bins
  st <- r$from + (seq_len(ncol(M)) - 1) * REP_WIN
  tabs[[r$name]] <- data.table(region = r$name, chrom = r$chrom, start = st, end = st + wlen,
                               masked = round(colSums(M), 4), t(round(M, 4)))
  w <- wlen / sum(wlen)
  summ[[r$name]] <- data.table(region = r$name, mb = round((r$to - r$from) / 1e6, 1),
                               masked_pct = round(100 * sum(colSums(M) * w), 1),
                               young_0_5_pct = round(100 * sum(M["0-5", ] * w), 1),
                               old_over20_pct = round(100 * sum((M["20-30", ] + M[">30", ]) * w), 1))
}

# all windows of all regions in one table, and the summary
if (length(tabs)) {
  out <- rbindlist(tabs, fill = TRUE); setnames(out, make.names(names(out), unique = TRUE))
  fwrite(out, file.path(OUTDIR, "repeats_by_window_all.tsv"), sep = "\t")
  cat("\n% of each region masked, and in the youngest (0-5%) and oldest (>20%) divergence bins:\n")
  print(rbindlist(summ))
}
cat("\nwrote", written, file.path(OUTDIR, "repeats_by_window_all.tsv"), sep = "\n  "); cat("\n")
