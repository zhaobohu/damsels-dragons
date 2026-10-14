# Statistics for Fig. 1B/1C: odonate vs mayfly displacement and ARI over the shared divergence range
# (console output only; optional argument overrides PAIRS)

library(data.table)

# paths - edit here
WORKDIR <- "/path/to/project/palaeoptera_synteny/pairwise_synteny"
PAIRS   <- "pairs_o2o.tsv"                     # from 02_synteny_pairs.R
K_MATCH <- 10          # odonate pairs matched to each mayfly pair by divergence

setwd(WORKDIR)
if (!is.na(commandArgs(trailingOnly = TRUE)[1]))
  PAIRS <- commandArgs(trailingOnly = TRUE)[1]
G_ODO <- "Odo-Odo"; G_EPH <- "Eph-Eph"

d <- fread(PAIRS)
need <- c("a", "b", "class", "dist", "ari", "nest", "na", "nb", "split90")
miss <- setdiff(need, names(d))
if (length(miss)) stop(basename(PAIRS), " is missing: ", paste(miss, collapse = ", "))
d[, displaced := 100 * (1 - nest)]
cat("[input] ", PAIRS, "   ", nrow(d), " pairs\n", sep = "")

# within-order pairs; shared range = at or below the deepest odonate pair
w_all <- d[class %chin% c(G_ODO, G_EPH)]
lim   <- d[class == G_ODO, max(dist)]
w     <- w_all[dist <= lim]
beyond <- w_all[class == G_EPH & dist > lim]

# 1. shared range and its divergence distribution
cat(sprintf("\n[1. shared range] dist <= %.4f (deepest odonate pair); %d mayfly pairs beyond it\n",
            lim, nrow(beyond)))
print(w[, .(pairs = .N, species = uniqueN(c(a, b)),
            d_q1 = round(quantile(dist, .25), 3), d_med = round(median(dist), 3),
            d_q3 = round(quantile(dist, .75), 3)), by = class][order(class)])

# 2. pooled medians, with and without the shared-range bound
summ <- function(x) x[, .(pairs = .N,
  disp_q1 = round(quantile(displaced, .25), 2), disp_med = round(median(displaced), 2),
  disp_q3 = round(quantile(displaced, .75), 2), disp_max = round(max(displaced), 2),
  ari_q1 = round(quantile(ari, .25), 3), ari_med = round(median(ari), 3),
  ari_q3 = round(quantile(ari, .75), 3)), by = class][order(class)]
cat("\n[2. pooled, shared range]\n"); print(summ(w))
cat(sprintf("  ratio of displacement medians (Eph/Odo): %.1f\n",
            w[class == G_EPH, median(displaced)] / w[class == G_ODO, median(displaced)]))
cat("\n[2b. all within-order pairs, i.e. WITHOUT the bound - for the text]\n")
print(summ(w_all))
cat(sprintf("  ratio without the bound: %.1f   (mayfly divergence median %.3f vs odonate %.3f)\n",
            w_all[class == G_EPH, median(displaced)] / w_all[class == G_ODO, median(displaced)],
            w_all[class == G_EPH, median(dist)], w_all[class == G_ODO, median(dist)]))

# 3. each mayfly pair against its K_MATCH nearest odonate pairs in divergence
# (gap = divergence distance to the furthest matched odonate pair)
o <- w[class == G_ODO]; e <- w[class == G_EPH]
mt <- list()
for (i in seq_len(nrow(e))) {
  nb <- o[order(abs(o$dist - e$dist[i]))][seq_len(min(K_MATCH, nrow(o)))]
  mt[[i]] <- data.table(a = e$a[i], b = e$b[i], dist = round(e$dist[i], 3),
                        eph_disp = round(e$displaced[i], 2),
                        odo_disp = round(median(nb$displaced), 2),
                        ratio    = round(e$displaced[i] / median(nb$displaced), 1),
                        disp_exceeds_all = e$displaced[i] > max(nb$displaced),
                        eph_ari  = round(e$ari[i], 3),
                        odo_ari  = round(median(nb$ari), 3),
                        ari_below_all = e$ari[i] < min(nb$ari),
                        gap = round(max(abs(nb$dist - e$dist[i])), 3))
}
mt <- rbindlist(mt)

cat(sprintf("\n[3. divergence-matched comparison, %d nearest odonate pairs each]\n", K_MATCH))
cat(sprintf("  displacement higher than EVERY matched odonate pair: %d of %d mayfly pairs\n",
            sum(mt$disp_exceeds_all), nrow(mt)))
cat(sprintf("  ARI lower than EVERY matched odonate pair:           %d of %d\n",
            sum(mt$ari_below_all), nrow(mt)))
cat(sprintf("  ratio of displacement to matched odonate median: %.1f (median), %.1f-%.1f\n",
            median(mt$ratio), min(mt$ratio), max(mt$ratio)))
cat(sprintf("  matching quality: furthest matched odonate pair is %.3f away at worst\n",
            max(mt$gap)))
ex <- mt[!(disp_exceeds_all)]
if (nrow(ex)) {
  cat("  mayfly pairs NOT exceeding every matched odonate pair:\n")
  print(ex[, .(a, b, dist, eph_disp, odo_disp, gap)])
  cat("  (shallow pairs where neither order has rearranged are expected here)\n")
}

# sensitivity of the count to K
cat("\n  sensitivity to K:\n")
for (kk in c(5, 10, 20, 40)) {
  if (kk > nrow(o)) next
  n_ex <- 0L
  for (i in seq_len(nrow(e))) {
    nb <- o[order(abs(o$dist - e$dist[i]))][seq_len(kk)]
    n_ex <- n_ex + (e$displaced[i] > max(nb$displaced))
  }
  cat(sprintf("    K = %2d: %d of %d\n", kk, n_ex, nrow(e)))
}

# 4. species level: median over each species' pairs
long <- rbind(w[, .(class, sp = a, displaced, ari)], w[, .(class, sp = b, displaced, ari)])
spm  <- long[, .(n_pairs = .N, disp = median(displaced), ari = median(ari)), by = .(class, sp)]
cat("\n[4. species medians]\n")
print(spm[order(class, -disp)][, .(class, sp, n_pairs, disp = round(disp, 2), ari = round(ari, 3))])
o <- spm[class == G_ODO]; e <- spm[class == G_EPH]
cat(sprintf("\n  displacement: highest odonate %.2f (%s)   lowest mayfly %.2f (%s)   %s\n",
            o[which.max(disp), disp], o[which.max(disp), sp],
            e[which.min(disp), disp], e[which.min(disp), sp],
            if (max(o$disp) < min(e$disp)) "NO OVERLAP" else "overlap"))
cat(sprintf("  ARI:          lowest odonate %.3f (%s)   highest mayfly %.3f (%s)   %s\n",
            o[which.min(ari), ari], o[which.min(ari), sp],
            e[which.max(ari), ari], e[which.max(ari), sp],
            if (min(o$ari) > max(e$ari)) "NO OVERLAP" else "overlap"))
cat("  ARI is reduced by chromosome-number change as well as by mixing, so a\n",
    "  fused or fissioned odonate (Tanypteryx) can sit below a conserved mayfly.\n",
    "  That is the reason both metrics are shown.\n", sep = "")

# 5. how localised the mixing is (split90)
cat("\n[5. fraction of chromosomes with < 90% of markers on one partner, favourable direction]\n")
print(w[, .(pairs = .N, split90_med = round(median(split90), 2),
            split90_min = round(min(split90), 2), split90_max = round(max(split90), 2)),
        by = class][order(class)])
cat("  low values with high displacement = mixing confined to a few chromosomes;\n",
    "  values near 1 = genome-wide.\n", sep = "")

# 6. Mann-Whitney on species medians (console only)
cat("\n[6. Mann-Whitney on species medians - CONSOLE ONLY, not for the manuscript]\n")
for (v in c("disp", "ari")) {
  t <- suppressWarnings(wilcox.test(o[[v]], e[[v]]))
  cat(sprintf("  %-4s W = %g   p = %.2g\n", v, t$statistic, t$p.value))
}
cat("  species within an order are phylogenetically non-independent; with complete\n",
    "  separation W is at its floor and p reflects n. Report the non-overlap.\n", sep = "")

# 7. caveats to carry into the text
cat("\n[7. caveats]\n",
    "  - both metrics read chromosome membership only: intrachromosomal\n",
    "    rearrangement is invisible\n",
    "  - a fission or fusion present in only one lineage of a pair is absorbed by\n",
    "    the favourable-direction rule; displacement is therefore a floor on\n",
    "    rearrangement in BOTH orders\n",
    "  - the relative substitution rate of the two orders is untested; a rate\n",
    "    difference would shift the comparison in an unknown direction\n", sep = "")
