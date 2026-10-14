# Prints every number in the neo-Y read-depth validation text, in manuscript order (reads genes.tsv, writes nothing)
library(data.table)

g <- fread("genes.tsv")

classes <- c("autosomal", "ancestral_X", "neoX_conserved", "neoX_lost")
labels  <- c("autosomal", "ancestral X", "retained (neoX_conserved)", "degenerated (neoX_lost)")
MIN_DEPTH <- 0.8    # flag: depth at or near autosomal level ...
HET_MULT  <- 1.5    # ... and heterozygosity above 1.5 x the autosomal median

fmt <- function(n) format(n, big.mark = ",")

# median normalized depth, all genes per class
cat("== Supp Results, paragraph 1: median normalized CDS depth (all genes) ==\n")
for (i in seq_along(classes)) {
  x <- g[class == classes[i]]
  cat(sprintf("  %-27s %.2f  (n = %s)\n", labels[i], median(x$norm_depth), fmt(nrow(x))))
}

# heterozygosity, genes with >= 300 callable CDS bases only (het_kb is NA below that)
h <- g[!is.na(het_kb)]
aut_het <- median(h[class == "autosomal", het_kb])
thr <- HET_MULT * aut_het

cat("\n== Supp Results, paragraph 1: heterozygosity (>= 300 callable CDS bases) ==\n")
cat(sprintf("  autosomal median het sites per callable kb: %.2f  (flag threshold %.2f)\n",
            aut_het, thr))
for (i in seq_along(classes)) {
  x <- h[class == classes[i]]
  cat(sprintf("  %-27s median %.2f; %.0f%% with no het site  (n = %s)\n",
              labels[i], median(x$het_kb), 100 * mean(x$het_n == 0), fmt(nrow(x))))
}

# flag: near-autosomal depth and high heterozygosity, applied to every class
h[, flag := norm_depth > MIN_DEPTH & het_kb > thr]

cat("\n== Flagged as candidate collapsed or unassembled gametologs ==\n")
for (i in seq_along(classes)) {
  x <- h[class == classes[i]]
  k <- sum(x$flag)
  cat(sprintf("  %-27s %d of %s (%.1f%%)\n", labels[i], k, fmt(nrow(x)), 100 * k / nrow(x)))
}

# degenerated vs ancestral X (the false-positive rate of the criterion)
deg <- h[class == "neoX_lost"]
anc <- h[class == "ancestral_X"]
tab <- matrix(c(sum(deg$flag), sum(!deg$flag), sum(anc$flag), sum(!anc$flag)), nrow = 2,
              dimnames = list(c("flagged", "not flagged"), c("degenerated", "ancestral X")))
ft <- fisher.test(tab)

cat("\n== Fisher's exact test, degenerated vs ancestral X ==\n")
print(tab)
cat(sprintf("  p = %.2g; odds ratio %.2f, 95%% CI %.2f-%.2f\n",
            ft$p.value, ft$estimate, ft$conf.int[1], ft$conf.int[2]))

# numbers for the manuscript placeholders
cat("\n== Fill-ins ==\n")
cat("Main Results: \"([deg %] against [anc %]; Fisher's exact test p = [p]; ...)\"\n")
cat(sprintf("  %.1f%% against %.1f%%; p = %.2g\n",
            100 * mean(deg$flag), 100 * mean(anc$flag), ft$p.value))
cat("Supp Results, paragraph 2:\n")
cat(sprintf("  %d of %s degenerated genes (%.1f%%) ... against %d of %s ancestral X genes (%.1f%%)\n",
            sum(deg$flag), fmt(nrow(deg)), 100 * mean(deg$flag),
            sum(anc$flag), fmt(nrow(anc)), 100 * mean(anc$flag)))
cat(sprintf("  Fisher's exact test p = %.2g; odds ratio %.2f, 95%% CI %.2f-%.2f\n",
            ft$p.value, ft$estimate, ft$conf.int[1], ft$conf.int[2]))
ret <- h[class == "neoX_conserved"]
cat(sprintf("  %d of the %s retained genes flagged\n", sum(ret$flag), fmt(nrow(ret))))
