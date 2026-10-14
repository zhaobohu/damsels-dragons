# Fig. S3: CDS depth by class (a), male-specific 31-mers by class (b), chr12 CDS depth by k (c); one svg each
library(data.table)

ROOT     <- "../.."                                   # repository root, prepended to every path below
RUNDEPTH <- "microchromosome/cds_depth"               # mchrom_k{3,2,1}_depth_per_gene.tsv
RUNKMER  <- "microchromosome/sex_kmers/kmer_sex_out"  # tables/windows_{num,den}*.tsv
OUTDIR   <- "."

# depth panels (a, c)
MIN_BP <- 300                 # minimum total CDS per gene
SRC_BP <- 300                 # the --min-cds-bp 03_mchrom_ratio.R was run with
KS     <- c(3, 2, 1)          # ladder rungs, drawn left to right
NBOOT  <- 2000

# k-mer panel (b): 100 kb windows pooled into rolling windows
ROLL_MB <- 2                  # rolling window width, Mb
STEP_MB <- 0.5                # rolling window step, Mb
BP_FILE <- "gene_pairs_ds/breakpoint.txt"             # fusion junction on the X

# compact letter displays (Games-Howell); panel b's use non-overlapping windows, panel c has none
PTHRESH    <- 0.05

set.seed(1)

W    <- 0.33      # half-width of the boxplots
SW   <- 0.42      # half-width of the beeswarm
STEP <- 0.013     # max spacing between adjacent points in a row
NMAX <- 1200      # points drawn per column; larger sets are subsampled

# sizes: CEX text, points and margins; LWD all lines (medians 2.2x); TCL tick length in text lines
CEX   <- 2
LWD   <- 2.5
TCL   <- -0.6
ALPHA <- 0.7      # point opacity

classes <- c("auto", "neoX_arm", "ancX_arm", "neoY", "micro")
labs    <- c("autosomal", "neo-X", "ancestral X", "neo-Y", "micro")
cols    <- c("#BBBBBB", "#66A61E", "#7570B3", "#B3717C", "#e6ab02")

MICRO <- "OZ369553.1"
XCHR  <- "OZ369554.1"
YCHR  <- "OZ369555.1"
AUTO  <- c("OZ369542.1","OZ369543.1","OZ369544.1","OZ369545.1","OZ369546.1","OZ369547.1",
           "OZ369548.1","OZ369549.1","OZ369550.1","OZ369551.1","OZ369552.1")

# paths under ROOT (empty parts skipped) and input check
p <- function(...) {
  parts <- c(path.expand(ROOT), ...)
  do.call(file.path, as.list(parts[nzchar(parts)]))
}
need <- function(f) {
  if (!file.exists(f))
    stop("not found: ", f,
         "\n  working directory: ", getwd(),
         "\n  a relative ROOT is resolved from there, not from the script's location",
         "\n  check ROOT, RUNDEPTH and RUNKMER at the top of this script")
  f
}

# fusion junction, required packages, CDS filter check
BP <- as.numeric(readLines(need(p(BP_FILE))))
for (pkg in c("rstatix", "multcompView"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(pkg, ' is not installed. Run install.packages("', pkg, '").')
if (MIN_BP < SRC_BP)
  warning("MIN_BP (", MIN_BP, ") is below the ", SRC_BP, " bp filter already applied by ",
          "03_mchrom_ratio.R; genes under ", SRC_BP, " bp are absent from the input and ",
          "cannot be recovered here.", call. = FALSE)

# per-gene depth for one k, genes with >= MIN_BP of CDS
load_k <- function(k) {
  d <- fread(need(p(RUNDEPTH, sprintf("mchrom_k%d_depth_per_gene.tsv", k))))
  d[cds_bp >= MIN_BP]
}
cat(sprintf("[depth] %d genes with a class and >= %d bp CDS\n", nrow(load_k(1)), MIN_BP))

# Games-Howell compact letters (group names must not contain "-")
cld_of <- function(values, groups, lv) {
  d <- data.frame(v = values, g = factor(as.character(groups), levels = lv))
  d <- d[is.finite(d$v), ]
  if (nlevels(droplevels(d$g)) < 2) return(setNames(rep("", length(lv)), lv))
  gh <- rstatix::games_howell_test(d, v ~ g)
  pv <- setNames(gh$p.adj, paste(gh$group1, gh$group2, sep = "-"))
  L  <- multcompView::multcompLetters(pv, threshold = PTHRESH)$Letters
  list(letters = L[lv], gh = gh,
       welch = rstatix::welch_anova_test(d, v ~ g))
}

# one panel: a beeswarm and a box per column, capped at ymax; rows above the data (top down: letters,
# n, bold `tops` labels) get headroom only when switched on
panel <- function(cols_v, tops, labels, colours, ylab, ticks, ymax,
                  unit = "genes", cld = NULL, show_n = TRUE, show_tops = F,
                  show_overflow = TRUE) {
  rows <- c(letters = !is.null(cld), n = show_n, tops = show_tops)
  nr   <- sum(rows)
  ypos <- setNames(rep(NA_real_, 3), c("letters", "n", "tops"))
  if (nr) ypos[names(rows)[rows]] <- ymax * (1.03 + 0.08 * rev(seq_len(nr) - 1))
  plot(NA, xlim = c(0.4, length(cols_v) + 0.6), ylim = c(0, ymax * (1.05 + 0.08 * nr)),
       axes = FALSE, xlab = "", ylab = ylab, yaxs = "i")
  axis(2, at = ticks, labels = ticks, las = 1, lwd = LWD)
  axis(1, at = seq_along(cols_v), labels = labels, tick = FALSE, line = 0.5, cex.axis = 0.85)
  box(bty = "l")
  nhi <- 0
  for (i in seq_along(cols_v)) {
    v <- cols_v[[i]][is.finite(cols_v[[i]])]
    if (!length(v)) next
    nhi <- nhi + sum(v > ymax)

    # beeswarm: at most NMAX points, none above the cap, each row of points spread sideways
    y <- if (length(v) > NMAX) sample(v, NMAX) else v
    y <- y[is.finite(y) & y <= ymax]
    bin <- round(y / (ymax / 90))
    xoff <- rep(0, length(y))
    for (bn in unique(bin)) {
      k <- which(bin == bn)
      k <- k[order(y[k])]
      step <- min(STEP, 2 * SW / length(k))       # per row, not per column
      xoff[k] <- (seq_along(k) - (length(k) + 1) / 2) * step
    }
    points(i + xoff, y, pch = 16,
           cex = if (length(y) < 200) 0.6 else 0.35,
           col = adjustcolor(colours[i], ALPHA))

    # box: quartiles, whiskers to the extreme values within 1.5 IQR, all stopped at the axis cap
    q <- quantile(v, c(0.25, 0.5, 0.75))
    iqr <- q[3] - q[1]
    b <- pmin(c(min(v[v >= q[1] - 1.5 * iqr]), q[1], q[2], q[3], max(v[v <= q[3] + 1.5 * iqr])), ymax)
    segments(i, b[1], i, b[2]); segments(i, b[4], i, b[5])
    if (b[1] < ymax) segments(i - W/2, b[1], i + W/2, b[1])
    if (b[5] < ymax) segments(i - W/2, b[5], i + W/2, b[5])
    rect(i - W, b[2], i + W, b[4], col = adjustcolor("white", 0.65), border = "black")
    segments(i - W, b[3], i + W, b[3], lwd = 2.2 * LWD, lend = 1)

    # labels above the column
    if (show_n)    text(i, ypos[["n"]], sprintf("n = %d", length(v)), cex = 0.72, col = "grey25")
    if (show_tops) text(i, ypos[["tops"]], tops[i], cex = 0.78, col = "grey15", font = 2)
    if (!is.null(cld) && nzchar(cld[i]))
      text(i, ypos[["letters"]], cld[i], font = 2, cex = 0.95)
  }
  # number of values above the cap (mtext sizes are absolute, hence par("cex"))
  if (show_overflow && nhi > 0) mtext(sprintf("%d %s above the axis", nhi, unit),
                     side = 3, line = -1.2, adj = 0.99, cex = 0.65 * par("cex"), col = "grey45")
}

# a: per-gene female/male CDS depth, each sex normalised to its autosomal median
A <- load_k(1)
A[, class := factor(class, levels = intersect(classes, unique(class)))]
am <- A[class == "auto", median(male)]; af <- A[class == "auto", median(fem)]
A[, r := (fem / af) / (male / am)]
lvA   <- levels(droplevels(A$class))

# a: per class, the gene ratios and the ratio of the class medians
vecsA <- list(); topsA <- c()
for (k in lvA) {
  vecsA[[length(vecsA) + 1]] <- A[class == k, r]
  topsA[k] <- (A[class == k, median(fem)] / af) / (A[class == k, median(male)] / am)
}
idxA  <- match(lvA, classes)

cat("\n[a: all CDS] class-level female/male (ratio of medians):\n")
print(data.table(class = lvA, n = sapply(vecsA, length), fem_over_male = round(topsA, 3)))

# a: Welch ANOVA, Games-Howell and letters
resA <- cld_of(A$r, A$class, lvA)
cldA <- resA$letters
cat("\n[a] Welch ANOVA:\n"); print(resA$welch)
cat("[a] Games-Howell:\n"); print(resA$gh[, c("group1","group2","estimate","p.adj")])
cat("[a] letters:", paste(lvA, cldA, sep = "="), "\n")

# b: male-specific (num) and single-copy (den) 31-mer counts per 100 kb window, with class
NW  <- 1e5
RW  <- round(ROLL_MB * 1e6 / NW)
RS  <- max(1, round(STEP_MB * 1e6 / NW))
sfx <- c("num1", "den1")     # single-copy track (unique in the assembly)
rdw <- function(f, v) fread(need(p(RUNKMER, "tables", f)), col.names = c("seq", "start", "end", v))
M <- merge(rdw(paste0("windows_", sfx[1], ".tsv"), "num"),
           rdw(paste0("windows_", sfx[2], ".tsv"), "den")[, .(seq, start, den)],
           by = c("seq", "start"))
M <- M[seq %in% c(AUTO, MICRO, XCHR, YCHR)][den > 0]
M[, class := fifelse(seq %in% AUTO, "auto",
             fifelse(seq == MICRO, "micro",
              fifelse(seq == YCHR, "neoY",
               fifelse(start < BP, "neoX_arm", "ancX_arm"))))]

# b: rolling windows within each chromosome and class, pooled (sum num / sum den, not a mean of
# percentages); overlapping windows are drawn (K), non-overlapping ones are tested (Kt)
K <- list(); Kt <- list()
for (d in split(M, by = c("seq", "class"))) {
  d <- d[order(start)]
  if (nrow(d) < RW) next
  pooled <- vapply(seq_len(nrow(d) - RW + 1), function(j)
    sum(d$num[j:(j + RW - 1)]) / sum(d$den[j:(j + RW - 1)]), numeric(1))
  i <- seq(1, nrow(d) - RW + 1, by = RS)
  K[[length(K) + 1]] <- data.table(class = d$class[i], v = pooled[i])
  i <- seq(1, nrow(d) - RW + 1, by = RW)
  Kt[[length(Kt) + 1]] <- data.table(class = d$class[i], v = pooled[i])
}
K <- rbindlist(K); Kt <- rbindlist(Kt)

# b: per class, the drawn windows and their median
lvK   <- intersect(classes, unique(K$class))
vecsK <- list(); medsK <- c()
for (k in lvK) {
  vecsK[[length(vecsK) + 1]] <- K[class == k, v]
  medsK <- c(medsK, median(K[class == k, v]))
}
idxK  <- match(lvK, classes)
amedK <- medsK[match("auto", lvK)]

cat("\n[b: male-specific k-mers] median % of the male's 31-mers the female lacks:\n")
print(data.table(class = lvK, n = sapply(vecsK, length), median_pct = round(medsK, 1),
                 fold_vs_auto = round(medsK / amedK, 2)))

# b: Welch ANOVA, Games-Howell and letters on the non-overlapping windows
resK <- cld_of(Kt$v, Kt$class, lvK)
cldK <- resK$letters
cat(sprintf("\n[b] letters computed on %d non-overlapping %g Mb windows (%d overlapping windows are drawn)\n",
            nrow(Kt), ROLL_MB, nrow(K)))
print(resK$welch)
cat("[b] Games-Howell:\n"); print(resK$gh[, c("group1","group2","estimate","p.adj")])
cat("[b] letters:", paste(lvK, cldK, sep = "="), "\n")

# c: for each k, chr12 shared genes normalised to autosomal shared genes; ratio of medians with a
# bootstrap CI over genes
vecsC <- list(); ladC <- list()
for (k in KS) {
  d  <- load_k(k)
  mi <- d[class == "micro" & shared == 1]
  au <- d[class == "auto"  & shared == 1]
  est <- (median(mi$fem) / median(au$fem)) / (median(mi$male) / median(au$male))
  bs <- replicate(NBOOT, {
    i <- sample.int(nrow(mi), replace = TRUE); j <- sample.int(nrow(au), replace = TRUE)
    (median(mi$fem[i]) / median(au$fem[j])) / (median(mi$male[i]) / median(au$male[j]))
  })
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
  vecsC[[length(vecsC) + 1]] <- (mi$fem / median(au$fem)) / (mi$male / median(au$male))
  ladC[[length(ladC) + 1]] <- data.table(k = k, n = nrow(mi), est = round(est, 3),
                                         lo = round(ci[1], 3), hi = round(ci[2], 3))
}

cat("\n[c: ladder] chr12 shared loci:\n")
print(rbindlist(ladC))
ymaxC <- max(1.5, ceiling(max(unlist(vecsC), na.rm = TRUE) * 10) / 10)

# draw the panels, one svg each (svglite keeps text editable; base svg() turns it into outlines)
written <- character()
for (pn in c("a.depth", "b.kmer", "c.sensitivity")) {
  f <- file.path(OUTDIR, paste0(pn, ".svg"))
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(f, width = 3000 / 190, height = 1900 / 190)
  } else {
    svg(f, width = 3000 / 190, height = 1900 / 190)
  }
  par(cex = CEX, lwd = LWD, tcl = TCL, mar = c(if (pn == "c.sensitivity") 5.2 else 4.5, 4.6, 2.6, 1),
      mgp = c(3.0, 0.7, 0))
  if (pn == "a.depth")
    panel(vecsA, sprintf("%.2f", topsA), labs[idxA], cols[idxA],
          "female / male CDS depth", c(0, 0.5, 1, 1.5, 2, 2.5, 3), 3,
          unit = "genes", cld = cldA,
          show_n = FALSE, show_overflow = FALSE)
  if (pn == "b.kmer")
    panel(vecsK, sprintf("%.0f%%", medsK), labs[idxK], cols[idxK],
          "proportion of male-specific 31-mers in 2Mb windows", seq(0, 1, 0.2), 1,
          unit = "windows", cld = cldK, show_n = FALSE)
  if (pn == "c.sensitivity") {
    # nothing above the data: the estimates and bootstrap CIs are printed above
    panel(vecsC, NULL,
          sprintf("%d", KS), rep("#e6ab02", length(KS)),
          "female / male CDS depth",
          seq(0, ymaxC, 0.5), ymaxC, unit = "genes", cld = NULL,
          show_n = FALSE, show_tops = FALSE, show_overflow = FALSE)
    mtext("A. juncea CDS shared with >=n other microchromosomes",
          side = 1, line = 3.6, cex = 0.8 * par("cex"))
  }
  invisible(dev.off())
  written <- c(written, f)
}

cat("\nwrote", written, sep = "\n  ")
cat("\n")
