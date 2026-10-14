# Supplementary comparison of the one-to-one (primary) and universal single-copy (robustness)
# marker sets from 02_synteny_pairs.R -> marker_set_comparison.tsv + agreement statistics

library(data.table)

# paths - edit here
WORKDIR <- "/path/to/project/palaeoptera_synteny/pairwise_synteny"
PRIM <- "pairs_o2o.tsv"                        # primary marker set
ROB  <- "pairs_universal.tsv"                  # robustness marker set
OUT  <- "marker_set_comparison.tsv"

setwd(WORKDIR)
G_ODO <- "Odo-Odo"; G_EPH <- "Eph-Eph"

# both pair tables are required
for (x in list(c(PRIM, "o2o"), c(ROB, "universal")))
  if (!file.exists(x[1]))
    stop(x[1], " does not exist. Create it with:\n  Rscript 02_synteny_pairs.R ", x[2],
         "\nThis script compares the two marker sets, so both are required.")

# pairs present in both tables; w = within-order pairs
p <- fread(PRIM); u <- fread(ROB)
cmp <- merge(p[, .(a, b, class, dist, n_o2o = n, na_o2o = na, nb_o2o = nb, ari_o2o = ari,
                   disp_o2o = 100 * (1 - nest))],
             u[, .(a, b, n_uni = n, na_uni = na, nb_uni = nb, ari_uni = ari,
                   disp_uni = 100 * (1 - nest))], by = c("a", "b"))
w <- cmp[class %chin% c(G_ODO, G_EPH)]

cat("[markers per pair] one-to-one median ", median(cmp$n_o2o),
    "   universal median ", median(cmp$n_uni), "\n", sep = "")

# Spearman correlation and median difference between the sets, per metric
cat("\n[agreement, within-order pairs]\n")
for (v in c("ari", "disp")) {
  x <- w[[paste0(v, "_o2o")]]; y <- w[[paste0(v, "_uni")]]
  r <- if (sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x, y, method = "spearman")
  cat(sprintf("  %-4s Spearman %s   median |diff| %.3f\n", v,
              if (is.na(r)) "n/a" else sprintf("%.3f", r), median(abs(x - y))))
}
cat("\n[chromosome counts] pairs where na or nb differs between sets: ",
    w[na_o2o != na_uni | nb_o2o != nb_uni, .N], " of ", nrow(w), "\n", sep = "")

# medians in the shared range (at or below the deepest odonate pair)
lim <- cmp[class == G_ODO, max(dist)]
s   <- w[dist <= lim]
tab <- s[, .(pairs = .N,
             disp_o2o = round(median(disp_o2o), 2), disp_uni = round(median(disp_uni), 2),
             ari_o2o = round(median(ari_o2o), 3),  ari_uni = round(median(ari_uni), 3)),
         by = class][order(class)]
cat(sprintf("\n[shared range, dist <= %.4f]\n", lim)); print(tab)
cat(sprintf("  displacement ratio Eph/Odo: one-to-one %.1f   universal %.1f\n",
            s[class == G_EPH, median(disp_o2o)] / s[class == G_ODO, median(disp_o2o)],
            s[class == G_EPH, median(disp_uni)] / s[class == G_ODO, median(disp_uni)]))

# species-level displacement (median over each species' pairs): do the orders overlap?
cat("\n[species-level displacement]\n")
for (col in c("disp_o2o", "disp_uni")) {
  l <- rbind(s[, .(class, sp = a, v = get(col))], s[, .(class, sp = b, v = get(col))])
  m <- l[, .(med = median(v)), by = .(class, sp)]
  cat(sprintf("  %-9s odonate max %.2f   mayfly min %.2f   %s\n", col,
              m[class == G_ODO, max(med)], m[class == G_EPH, min(med)],
              if (m[class == G_ODO, max(med)] < m[class == G_EPH, min(med)]) "no overlap" else "OVERLAP"))
}

# write the medians table (shared range and all within-order pairs)
out <- rbind(
  cbind(scope = "shared range", tab),
  cbind(scope = "all within-order",
        w[, .(pairs = .N,
              disp_o2o = round(median(disp_o2o), 2), disp_uni = round(median(disp_uni), 2),
              ari_o2o = round(median(ari_o2o), 3),  ari_uni = round(median(ari_uni), 3)),
          by = class][order(class)]))
fwrite(out, OUT, sep = "\t")
cat("\n[wrote] ", OUT, "\n", sep = "")

# the pairs whose displacement changes most between the sets
cat("\n[largest displacement changes]\n")
print(w[order(-abs(disp_o2o - disp_uni))][1:8, .(a, b, class, n_uni, n_o2o,
      disp_uni = round(disp_uni, 2), disp_o2o = round(disp_o2o, 2))])
