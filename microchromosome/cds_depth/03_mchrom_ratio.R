#!/usr/bin/env Rscript
# Gene-level female:male CDS depth by class (all CDS, and shared orthologs only), each sex normalised
# to its own autosomal median. Run where 02_depth_both_sexes.sh wrote its output:
#   Rscript 03_mchrom_ratio.R --male male.cds.q20.regions.bed.gz --female female.cds.q20.regions.bed.gz \
#     --classes gene_classes_mchrom.tsv [--min-cds-bp 300] [--nboot 2000] [--out mchrom]

suppressPackageStartupMessages(library(data.table))

# command-line arguments
a <- commandArgs(TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, a)
  if (is.na(i)) return(default)
  if (i == length(a)) stop("missing value for ", flag)
  a[i + 1L]
}
male_f  <- getarg("--male")
fem_f   <- getarg("--female")
cls_f   <- getarg("--classes")
min_bp  <- as.numeric(getarg("--min-cds-bp", "300"))
nboot   <- as.integer(getarg("--nboot", "2000"))
outpref <- getarg("--out", "mchrom")
if (any(vapply(list(male_f, fem_f, cls_f), is.null, TRUE)))
  stop("--male, --female and --classes are required")
for (f in c(male_f, fem_f, cls_f)) if (!file.exists(f)) stop("missing: ", f)

CLS <- c("auto", "neoX_arm", "ancX_arm", "neoX", "neoY", "micro")

# per-gene depth: length-weighted mean over that gene's CDS
read_depth <- function(path, label) {
  d <- fread(cmd = paste("gzip -cd", shQuote(path)), header = FALSE,
             col.names = c("chrom", "start", "end", "gene", "mean_depth"))
  d[, len := end - start]
  out <- d[, .(depth = sum(mean_depth * len) / sum(len), cds_bp = sum(len)),
           by = gene]
  setnames(out, "depth", label)
  out
}

# both sexes, genes with enough CDS, joined to the gene classes
M  <- read_depth(male_f, "male")
F_ <- read_depth(fem_f,  "fem")
D  <- merge(M, F_[, .(gene, fem)], by = "gene")
cat(sprintf("[depth] genes measured in both sexes: %d\n", nrow(D)))
D <- D[cds_bp >= min_bp]
cat(sprintf("[depth] after >= %.0f bp CDS: %d\n", min_bp, nrow(D)))

G <- fread(cls_f)
D <- merge(D, G[, .(gene, chrom, class, shared)], by = "gene")
cat(sprintf("[depth] genes with a class: %d\n", nrow(D)))
if (!nrow(D))
  stop("no gene IDs in common between the depth files and ", cls_f,
       " -- different namespaces?")
D[, class := factor(class, levels = intersect(CLS, unique(class)))]

# per-class medians, normalised to the autosomal medians of the same subset
by_class <- function(dt, tag) {
  auto_m <- dt[class == "auto", median(male, na.rm = TRUE)]
  auto_f <- dt[class == "auto", median(fem,  na.rm = TRUE)]
  if (!is.finite(auto_m) || !is.finite(auto_f) || auto_m == 0 || auto_f == 0)
    stop("no usable autosomal genes to normalise on, subset: ", tag)
  out <- dt[, .(n = .N,
                male_med = median(male, na.rm = TRUE),
                fem_med  = median(fem,  na.rm = TRUE)), by = class]
  out[, `:=`(subset = tag,
             male_norm = male_med / auto_m,
             fem_norm  = fem_med  / auto_f)]
  out[, fem_over_male := fem_norm / male_norm]
  out[order(match(as.character(class), CLS))]
}

# all CDS and shared orthologs
S <- D[shared == 1]
tab <- rbind(by_class(D, "all_cds"), by_class(S, "shared"))
setcolorder(tab, c("subset", "class", "n", "male_med", "fem_med",
                   "male_norm", "fem_norm", "fem_over_male"))

# bootstrap CI of the shared micro fem_over_male (resampling micro and auto genes)
# sample.int rather than sample(x): sample(x) on a length-1 vector draws from 1:x
i_mi <- which(S$class == "micro"); i_au <- which(S$class == "auto")
ci <- c(NA_real_, NA_real_)
if (length(i_mi) >= 5 && length(i_au) >= 5) {
  set.seed(1)
  bs <- replicate(nboot, {
    im <- i_mi[sample.int(length(i_mi), replace = TRUE)]
    ia <- i_au[sample.int(length(i_au), replace = TRUE)]
    (median(S$fem[im]) / median(S$fem[ia])) / (median(S$male[im]) / median(S$male[ia]))
  })
  ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
} else {
  warning("fewer than 5 shared micro or auto genes -- no bootstrap; the point ",
          "estimate should not be quoted without saying so")
}

# console summary and tables
cat("\n== CDS depth / autosomal median, by class (fem_over_male) ==\n")
print(tab)
pt <- tab[subset == "shared" & class == "micro", fem_over_male]
if (length(pt) && is.finite(ci[1]))
  cat(sprintf("\nmicro, shared orthologs: fem_over_male = %.3f  95%% CI %.3f - %.3f  (n = %d)\n",
              pt, ci[1], ci[2], length(i_mi)))
cat("\n1.0 = equal in both sexes (repeat/mismapping artifact)\n",
    "2.0 = male-hemizygous, X-like;  0.5 = female-hemizygous, ZW-like\n",
    "compare the micro rows against tables/depth_by_class.tsv (100-kb windows).\n", sep = "")

fwrite(tab, paste0(outpref, "_depth_by_class.tsv"), sep = "\t")
fwrite(D[, .(gene, chrom, class, shared, cds_bp, male, fem)],
       paste0(outpref, "_depth_per_gene.tsv"), sep = "\t")

# barplot of fem_over_male per class, one panel per subset
pdf(paste0(outpref, "_depth.pdf"), width = 10, height = 4.5)
op <- par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
for (tg in c("all_cds", "shared")) {
  dt <- if (tg == "shared") S else D
  am <- dt[class == "auto", median(male)]; af <- dt[class == "auto", median(fem)]
  lv <- levels(droplevels(dt$class))
  r <- c()
  for (k in lv) r[k] <- (median(dt$fem[dt$class == k]) / af) / (median(dt$male[dt$class == k]) / am)
  r[!is.finite(r)] <- NA                     # a zero male median would give Inf
  top <- max(2.2, max(r, na.rm = TRUE) * 1.1)
  bp <- barplot(r, ylim = c(0, top), las = 2,
                ylab = "fem / male (autosome-normalised)", main = tg)
  abline(h = c(0.5, 1, 2), lty = c(3, 2, 3), col = c("red", "grey40", "blue"))
  text(bp, r, ifelse(is.na(r), "", sprintf("%.2f", r)), pos = 3, cex = 0.8, xpd = NA)
}
par(op)
invisible(dev.off())

cat(sprintf("\nwrote %s_depth_by_class.tsv, %s_depth_per_gene.tsv, %s_depth.pdf\n",
            outpref, outpref, outpref))
