library(data.table)
library(beeswarm)
library(rstatix)
library(multcompView)

SMOOTH_WIN   <- 2000             # rolling-mean window in 5 kb bins -> 10 Mb
SMOOTH_STEP  <- 30               # step in 5 kb bins -> 150 kb
MALE_MIN     <- 0.15             # male-depth floor (x autosomal median)
STAT_BIN     <- 3e6              # bin size for the 5-category stats

# 5 kb window depths: mosdepth -Q 50 --by windows.map.bed on align.uniq.bam
# (read_depth/*/2.Aeshna_juncea_depth.sbatch)
fem <- fread("../../read_depth/female_SRR36191036/SRR36191036.uniq.regions.bed.gz")
mal <- fread("../../read_depth/male_ERR12340099/ERR12340099.uniq.regions.bed.gz")
setnames(fem, c("query","wstart","wend","mapfrac","depth"))
setnames(mal, c("query","wstart","wend","mapfrac","depth"))

# chromosome names
chroms <- c(
  "OZ369542.1"="1","OZ369543.1"="2","OZ369544.1"="3","OZ369545.1"="4",
  "OZ369546.1"="5","OZ369547.1"="6","OZ369548.1"="7","OZ369549.1"="8",
  "OZ369550.1"="9","OZ369551.1"="10","OZ369552.1"="11","OZ369553.1"="12",
  "OZ369554.1"="X","OZ369555.1"="Y"
)
lev <- c(as.character(1:12), "X", "Y")
mapchr <- function(d) {
  d <- d[query %in% names(chroms)]
  d[, query := factor(chroms[query], levels = lev)][]
}
fem <- mapchr(fem); mal <- mapchr(mal)

# merge sexes and normalize each by its autosomal median depth
d  <- merge(fem, mal, by = c("query","wstart","wend"), suffixes = c(".f",".m"))
mf <- d[query %in% as.character(1:12), median(depth.f)]
mm <- d[query %in% as.character(1:12), median(depth.m)]
d[, `:=`(depth.f = depth.f / mf, depth.m = depth.m / mm)]

# female:male depth ratio, scaled by the mean autosomal ratio
d <- d[depth.m >= MALE_MIN]
d[, ratio := depth.f / depth.m]
ra <- d[query %in% as.character(1:12), mean(ratio)]
d[, ratio := ratio / ra]
setorder(d, query, wstart)

# rolling mean per chromosome; chromosomes with fewer than SMOOTH_WIN bins keep raw ratios
out <- NULL
for (ch in levels(d$query)) {
  x <- d[query == ch]
  n <- nrow(x)
  if (n < 1) next
  if (n < SMOOTH_WIN) {
    out <- rbind(out, data.frame(chr = x$query, mid = (x$wstart + x$wend) / 2,
                                 mean = x$ratio))
    next
  }
  i <- seq(1, n - SMOOTH_WIN + 1, by = SMOOTH_STEP)
  out <- rbind(out, data.frame(
    chr  = x$query[i],
    mid  = (x$wstart[i] + x$wend[i + SMOOTH_WIN - 1]) / 2,
    mean = sapply(i, function(j) mean(x$ratio[j:(j + SMOOTH_WIN - 1)]))
  ))
}

# split the X at the fusion breakpoint
BP_FILE <- "../../gene_pairs_ds/breakpoint.txt"
X_CUT <- as.numeric(readLines(BP_FILE))
lev2 <- c(as.character(1:12), "neo-X", "anc-X", "Y")
out$chr <- as.character(out$chr)
out$chr[out$chr == "X" & out$mid <  X_CUT] <- "neo-X"
out$chr[out$chr == "X" & out$mid >= X_CUT] <- "anc-X"
out$chr <- droplevels(factor(out$chr, levels = lev2))

# plot
smooth_mb <- SMOOTH_WIN  * 5e3 / 1e6
step_kb <- SMOOTH_STEP * 5e3 / 1e3
plot_title <- sprintf("female:male depth in %g Mb windows (%g kb step)",
                       smooth_mb, step_kb)

svg("depth.ratio.boxplot.svg", width = 22, height = 7)
lv <- levels(out$chr); n <- length(lv)
bg <- c("white","grey93")[seq_along(lv) %% 2 + 1]
bg[lv == "neo-X"] <- adjustcolor("#7570B3", 0.4)
bg[lv == "anc-X"] <- adjustcolor("#66A61E", 0.4)
bg[lv == "Y"]     <- adjustcolor("#B3717C", 0.4)
plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0, 3.5), xaxs = "i", yaxs = "i",
     xaxt = "n", xlab = "chromosome", ylab = "female:male depth",
     mgp = c(3, 3.5, 0), tcl = -2,
     main = plot_title)
rect(seq_len(n) - 0.5, par("usr")[3], seq_len(n) + 0.5, par("usr")[4],
     col = bg, border = NA)
axis(1, at = seq_len(n), labels = lv, tick = FALSE)
beeswarm(out$mean ~ out$chr, add = TRUE, cex = 0.18, pch = 16,
         spacing = 1.4, col = adjustcolor("black", 0.9))
box()
dev.off()

# mean ratio in non-overlapping STAT_BIN bins, grouped into 5 categories
micro <- "12"
b <- d[, .(ratio = mean(ratio), mid = mean((wstart + wend)/2)),
       by = .(query, bin = wstart %/% STAT_BIN)]
cat5 <- function(chr, mid) {
  chr <- as.character(chr)
  g <- ifelse(chr %in% as.character(1:12), "auto", chr)
  g[chr %in% micro] <- "micro"
  g[chr == "X" & mid <  X_CUT] <- "neoX"
  g[chr == "X" & mid >= X_CUT] <- "ancX"
  g[chr == "Y"]                <- "Y"
  factor(g, c("auto","micro","neoX","ancX","Y"))
}
b <- b[!is.na(cat5(query, mid))][, grp := cat5(query, mid)]

# pairwise Games-Howell tests
gl <- levels(b$grp)
gh <- games_howell_test(b, ratio ~ grp)                   # means (FWER)

# Cohen's d with bootstrap 95% CIs, all pairs
set.seed(1)
cd <- cohens_d(b, ratio ~ grp, ci = TRUE, ci.type = "perc", nboot = 2000)
cat("\n[Cohen's d, all pairwise (bootstrap 95% CIs)]:\n")
print(cd[, c("group1", "group2", "effsize", "conf.low", "conf.high", "magnitude")],
      row.names = FALSE)

# compact letters from Games-Howell
pv  <- setNames(gh$p.adj, paste(gh$group1, gh$group2, sep = "-"))
cld <- multcompLetters(pv, threshold = 0.05)$Letters

# 5-category boxplot with letters
disp <- c(auto="auto", micro="micro", neoX="neo-X", ancX="anc-X", Y="Y")
b[, grp2 := factor(disp[as.character(grp)], disp)]
lvb <- levels(b$grp2); k <- length(lvb)
bg <- c("white","grey93")[seq_len(k) %% 2 + 1]
bg[lvb=="neo-X"] <- "#fddbc7"; bg[lvb=="anc-X"] <- "#f4a582"; bg[lvb=="Y"] <- "#92c5de"
plot(NA, xlim = c(0.5, k+0.5), ylim = c(0, 3), xaxs="i", yaxs="i", xaxt="n",
     xlab = "", ylab = "female:male depth", main = "")
rect(seq_len(k)-0.5, par("usr")[3], seq_len(k)+0.5, par("usr")[4], col=bg, border=NA)
axis(1, at = seq_len(k), labels = lvb, tick = FALSE)
beeswarm(ratio ~ grp2, data = b, add = TRUE, cex = 0.2, pch = 16,
         col = adjustcolor("black", 0.25))
boxplot(ratio ~ grp2, data = b, add = TRUE, outline = FALSE, col = NA,
        lwd = 2, border = "black", xaxt = "n", yaxt = "n")
tt   <- b[, .(top = max(ratio)), by = grp2]
tops <- tt$top[match(lvb, as.character(tt$grp2))]
text(seq_len(k), pmin(tops + 0.15, 2.9), cld[names(disp)], font = 2)
box()

# omnibus Welch ANOVA (printed: F, df, p)
welch_anova_test(b, ratio ~ grp)

# forest plot of Games-Howell estimates + CIs
g   <- gh[order(gh$estimate), ]
lab <- paste(g$group1, g$group2, sep = " vs ")
plot(g$estimate, seq_along(lab), pch = 16,
     xlim = range(c(g$conf.low, g$conf.high, 0)), yaxt = "n",
     ylab = "", xlab = expression(Delta~"female:male depth (group2 - group1)"))
segments(g$conf.low, seq_along(lab), g$conf.high, seq_along(lab), lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
axis(2, at = seq_along(lab), labels = lab, las = 1)
