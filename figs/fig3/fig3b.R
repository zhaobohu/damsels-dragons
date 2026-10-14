library(data.table)
library(beeswarm)
library(rstatix)
library(multcompView)

# read the RepeatMasker output
dat <- as.data.frame(fread(
  "../../annotation/repeats/Aeshna_juncea/Aeshna.juncea.clean.fna.out.gz",
  header = FALSE, skip = 3, fill = TRUE, comment.char = "", quote = ""
))
colnames(dat) <- c(
  "SW.score","perc.div","perc.del","perc.ins","query","qstart","qend",
  "left","strand","repeat","class.family","sstart","send","left2","id","flag"
)

# drop hits flagged "*" (overlapped by a higher-scoring hit) and hits scoring < 225,
# except simple repeats / low complexity
dat <- dat[dat$flag != "*", ]
dat <- dat[dat$SW.score >= 225 |
             dat$class.family %in% c("Simple_repeat","Low_complexity"), ]

# chromosome names
chroms <- c(
  "OZ369542.1"="1","OZ369543.1"="2","OZ369544.1"="3","OZ369545.1"="4",
  "OZ369546.1"="5","OZ369547.1"="6","OZ369548.1"="7","OZ369549.1"="8",
  "OZ369550.1"="9","OZ369551.1"="10","OZ369552.1"="11","OZ369553.1"="12",
  "OZ369554.1"="X","OZ369555.1"="Y"
)
chroms <- factor(chroms, levels = chroms)
dat <- dat[dat$query %in% names(chroms), ]
dat$query <- chroms[dat$query]
dat$query <- factor(dat$query, levels(chroms))

# make global coordinates
sizes <- as.data.frame(fread("../../data/genomes/Aeshna.juncea.clean.masked.fna.fai"))
sizes <- sizes[sizes$V1 %in% names(chroms), ]
sizes$V1 <- chroms[sizes$V1]
len <- setNames(sizes$V2, sizes$V1)
off <- c(0, head(cumsum(len), -1))
names(off) <- levels(dat$query)
o <- off[as.character(dat$query)]
dat$gstart <- dat$qstart + o
dat$gend <- dat$qend + o
dat <- dat[, c("query","gstart","gend")]

# 5 kb windows along each chromosome (global coordinates)
w <- 5000
wins <- list()
for (chr in levels(dat$query)) {
  s <- seq(0, len[chr] - 1, by = w) + off[chr]
  wins[[chr]] <- data.table(query = chr, wstart = s, wend = pmin(s + w, off[chr] + len[chr]))
}
wins <- rbindlist(wins)

# repeat bp and repeat proportion per window
rp <- as.data.table(dat)
setkey(wins, query, wstart, wend)
ov <- foverlaps(rp, wins, by.x = c("query","gstart","gend"), nomatch = NULL)
ov[, olen := pmin(gend, wend) - pmax(gstart, wstart)]
agg <- ov[, .(repeat.bp = sum(olen)), by = .(query, wstart, wend)]
wins <- merge(wins, agg, by = c("query","wstart","wend"), all.x = TRUE)
wins[is.na(repeat.bp), repeat.bp := 0]
wins[, `:=`(prop = repeat.bp / (wend - wstart), mid = (wstart + wend)/2)]
wins[, query := factor(query, levels(dat$query))]

# rolling mean per chromosome; chromosomes with fewer than SMOOTH_WIN windows are dropped
SMOOTH_WIN  <- 2000   # windows -> 10 Mb
SMOOTH_STEP <- 30     # windows -> 150 kb
out <- NULL
for (ch in levels(wins$query)) {
  x <- wins[query == ch]
  if (nrow(x) < SMOOTH_WIN) next
  i <- seq(1, nrow(x) - SMOOTH_WIN + 1, by = SMOOTH_STEP)
  out <- rbind(out, data.frame(
    chr  = x$query[i],
    mid  = (x$wstart[i] + x$wend[i + SMOOTH_WIN - 1]) / 2,
    mean = sapply(i, function(j) mean(x$prop[j:(j + SMOOTH_WIN - 1)]))
  ))
}

# split the X at the inferred breakpoint (BP is a local X coordinate)
BP <- as.numeric(readLines("../../gene_pairs_ds/breakpoint.txt"))
lev2 <- c(as.character(1:12), "neo-X", "anc-X", "Y")
loc <- out$mid - off[["X"]]
out$chr <- as.character(out$chr)
out$chr[out$chr == "X" & loc <  BP] <- "neo-X"
out$chr[out$chr == "X" & loc >= BP] <- "anc-X"
out$chr <- factor(out$chr, levels = lev2)

# plot
svg("repeat.boxplot.svg", width = 22, height = 7)
lev <- levels(out$chr); n <- length(lev)
bg <- c("white","grey93")[seq_along(lev) %% 2 + 1]
bg[lev == "neo-X"] <- adjustcolor("#7570B3", 0.4)
bg[lev == "anc-X"] <- adjustcolor("#66A61E", 0.4)
bg[lev == "Y"]     <- adjustcolor("#B3717C", 0.4)
plot(NA, xlim = c(0.5, n + 0.5), ylim = c(0.1, 0.9), xaxs = "i", yaxs = "i",
     xaxt = "n", xlab = "chromosome", ylab = "repeat proportion",
     mgp = c(3, 3.5, 0), tcl = -2,
     main = "repeat proportion in 10 Mb windows (150 kb step)")
rect(seq_len(n) - 0.5, par("usr")[3], seq_len(n) + 0.5, par("usr")[4],
     col = bg, border = NA)
axis(1, at = seq_len(n), labels = lev, tick = FALSE)
beeswarm(out$mean ~ out$chr, add = TRUE, cex = 0.3, pch = 16,
         col = adjustcolor("black", 0.5))
box()
dev.off()

# mean proportion in non-overlapping 3 Mb bins, grouped into 5 categories
B <- 3e6 # bin size for stats
micro <- "12" # microchromosome
b <- wins[, .(prop = mean(prop), mid = mean(mid)),
          by = .(query, bin = (wstart - off[as.character(query)]) %/% B)]
cat5 <- function(chr, mid) {
  chr <- as.character(chr)
  loc <- mid - off[["X"]] # wins mid is GLOBAL -> local X
  g <- ifelse(chr %in% as.character(1:12), "auto", chr)
  g[chr %in% micro] <- "micro"
  g[chr == "X" & loc < BP] <- "neoX"
  g[chr == "X" & loc >= BP] <- "ancX"
  g[chr == "Y"] <- "Y"
  factor(g, c("auto","micro","neoX","ancX","Y"))
}
b <- b[!is.na(cat5(query, mid))][, grp := cat5(query, mid)]

# omnibus Welch ANOVA (printed: F, df, p) and pairwise Games-Howell tests
welch_anova_test(b, prop ~ grp)
gh <- games_howell_test(b, prop ~ grp) # pairwise (fwer)

# Cohen's d with bootstrap 95% CIs, all pairs
set.seed(1)
cd <- cohens_d(b, prop ~ grp, ci = TRUE, ci.type = "perc", nboot = 2000)
cat("\n[Cohen's d, all pairwise (bootstrap 95% CIs)]:\n")
print(cd[, c("group1", "group2", "effsize", "conf.low", "conf.high", "magnitude")],
      row.names = FALSE)

# compact letters from Games-Howell
pv  <- setNames(gh$p.adj, paste(gh$group1, gh$group2, sep = "-"))
cld <- multcompLetters(pv, threshold = 0.05)$Letters

# 5-category boxplot with letters
disp <- c(auto="auto", micro="micro", neoX="neo-X", ancX="anc-X", Y="Y")
b[, grp2 := factor(disp[as.character(grp)], disp)]
lv <- levels(b$grp2); k <- length(lv)
bg <- c("white","grey93")[seq_len(k) %% 2 + 1]
bg[lv=="neo-X"] <- "#fddbc7"; bg[lv=="anc-X"] <- "#f4a582"; bg[lv=="Y"] <- "#92c5de"
ymax <- max(b$prop) * 1.05
plot(NA, xlim = c(0.5, k+0.5), ylim = c(0, ymax), xaxs="i", yaxs="i", xaxt="n",
     xlab = "", ylab = "repeat proportion", main = "")
rect(seq_len(k)-0.5, par("usr")[3], seq_len(k)+0.5, par("usr")[4], col=bg, border=NA)
axis(1, at = seq_len(k), labels = lv, tick = FALSE)
beeswarm(prop ~ grp2, data = b, add = TRUE, cex = 0.4, pch = 16,
         col = adjustcolor("black", 0.25))
boxplot(prop ~ grp2, data = b, add = TRUE, outline = FALSE, col = NA,
        lwd = 2, border = "black", xaxt = "n", yaxt = "n")
tt   <- b[, .(top = max(prop)), by = grp2]
tops <- tt$top[match(lv, as.character(tt$grp2))]
text(seq_len(k), pmin(tops + 0.03*ymax, ymax*0.99), cld[names(disp)], font = 2)
box()

# forest plot of Games-Howell estimates + CIs
g   <- gh[order(gh$estimate), ]
lab <- paste(g$group1, g$group2, sep = " vs ")
plot(g$estimate, seq_along(lab), pch = 16,
     xlim = range(c(g$conf.low, g$conf.high, 0)), yaxt = "n",
     ylab = "", xlab = expression(Delta~"repeat proportion (group2 - group1)"))
segments(g$conf.low, seq_along(lab), g$conf.high, seq_along(lab), lwd = 2)
abline(v = 0, lty = 2, col = "grey50")
axis(2, at = seq_along(lab), labels = lab, las = 1)
