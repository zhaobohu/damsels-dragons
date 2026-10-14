# Out-of-X gene traffic (Supplementary Results), on orthogroups single-copy in all ten genomes
suppressMessages(library(data.table))

PG   <- "../../odonate_synteny/genespace/pangenes/Aeshna_juncea_pangenes.txt.gz"
AJ   <- "Aeshna_juncea";        AJ_X <- "OZ369554.1"; AJ_Y <- "OZ369555.1"
ANAX <- "Anax_junius";          ANAX_NEOX <- "OZ349811.1"; ANAX_ANCX <- "OZ349813.1"
SOMA <- "Somatochlora_arctica"; SOMA_NEOX <- "OZ407514.1"; SOMA_ANCX <- "OZ407525.1"
CONF <- c("PASS", "array", "NSOrtho"); OVERMERGE_CHR <- 3   # unlike fig4.R, NSOrtho counts as confident
#   revert to c("PASS", "array") to reproduce the syntenic-only version

# fusion breakpoint written by gene_pairs_ds/pairs.sbatch (103 Mb if the file is missing)
BP_FILE <- "../../gene_pairs_ds/breakpoint.txt"
BRK <- if (file.exists(BP_FILE)) {
  as.numeric(readLines(BP_FILE))
} else {
  warning("breakpoint.txt not found; falling back to 103e6")
  103e6
}
cat(sprintf("[breakpoint] using %.3f Mb (from %s)\n", BRK / 1e6, BP_FILE))

gid <- function(x) sub("\\.t[0-9]+$", "", x)

pg <- fread(PG); pg[, g := gid(id)]; pg[, start := as.numeric(start)]; pg[is.na(og), og := -pgID]

# juncea arm of each gene, and confident juncea copies per og per arm (wide table)
jloc <- unique(pg[genome == AJ, .(g, chr, start)])
jloc[, arm := fifelse(chr == AJ_X & start <  BRK, "neoX",
              fifelse(chr == AJ_X & start >= BRK, "ancX",
              fifelse(chr == AJ_Y,                "Y",
              fifelse(grepl("^OZ", chr),          "auto", "unplaced"))))]
jloc <- unique(jloc[, .(jg = g, arm, chr)])
jc   <- merge(pg[genome == AJ & flag %chin% CONF, .(g, og)], jloc, by.x = "g", by.y = "jg")
jcnt <- dcast(jc, og ~ arm, value.var = "g", fun.aggregate = uniqueN, fill = 0L)
for (cc in c("neoX","ancX","Y","auto","unplaced"))
  if (!cc %in% names(jcnt)) jcnt[, (cc) := 0L]

# over-merge filter: Anax chromosomes per og
og_nchr <- pg[genome == ANAX & grepl("^OZ", chr) & flag %chin% CONF,
              .(nchr = uniqueN(chr)), by = og]

# strict single-copy orthogroups: one confident gene in each of the ten genomes (the og spans all ten)
sc <- pg[flag %chin% CONF, .(n = uniqueN(g)), by = .(og, genome)][
         , .(ngen = .N, mx = max(n)), by = og][ngen == 10L & mx == 1L, og]
cat(sprintf("[strict set] %d single-copy orthogroups present in all 10 genomes\n",
            length(sc)))

# Soma OZ407514.1 split per gene by its og's Anax chromosome (neo-X homolog vs auto-S), as in fig4.R
og_achr <- pg[genome == ANAX & grepl("^OZ", chr) & flag %chin% CONF,
              .(N = uniqueN(g)), by = .(og, chr)][
              order(-N), .(achr = chr[1]), by = og]
sf <- merge(pg[genome == SOMA & chr == SOMA_NEOX, .(g, og)], og_achr, by = "og", all.x = TRUE)
sf[, blk := fifelse(!is.na(achr) & achr == ANAX_NEOX, "neoXhom", "autoS")]
soma_neox_genes <- sf[blk == "neoXhom", g]
soma_ancx_genes <- pg[genome == SOMA & chr == SOMA_ANCX, g]
cat(sprintf("[Soma fusion split] OZ407514.1: neo-X-homolog %d, auto-S excluded %d\n",
            length(soma_neox_genes), nrow(sf) - length(soma_neox_genes)))

# neo-X- and anc-X-source ogs: Anax ogs with a confident copy on the source chromosome
build_x <- function(anax_src, soma_genes) {
  focus <- pg[genome == ANAX & chr == anax_src & flag %chin% CONF, .(a = uniqueN(g)), by = og]
  scnt  <- pg[genome == SOMA & g %chin% soma_genes & flag %chin% CONF, .(s = uniqueN(g)), by = og]
  G <- focus
  for (x in list(scnt, jcnt, og_nchr)) G <- merge(G, x, by = "og", all.x = TRUE)
  for (cc in c("s","neoX","ancX","Y","auto","unplaced","nchr"))
    G[is.na(get(cc)), (cc) := 0L]
  G[, keep := nchr <= OVERMERGE_CHR & og %chin% sc][]
}
G_neoX <- build_x(ANAX_NEOX, soma_neox_genes)
G_ancX <- build_x(ANAX_ANCX, soma_ancx_genes)

# autosomal-source ogs: Anax ogs on any placed autosome
anax_chrs  <- pg[genome == ANAX & grepl("^OZ", chr), unique(chr)]
anax_autos <- setdiff(anax_chrs, c(ANAX_NEOX, ANAX_ANCX))

# Anax focal chromosome and copy count per og (largest wins; over-merged ogs are dropped later)
af <- pg[genome == ANAX & chr %chin% anax_autos & flag %chin% CONF,
         .(N = uniqueN(g)), by = .(og, achr = chr)][order(og, -N)]
focus_auto <- af[!duplicated(og), .(og, achr, a = N)]

# expected home: the modal juncea and Soma chromosome among the ogs of each Anax autosome
ja <- merge(pg[genome == AJ   & flag %chin% CONF & grepl("^OZ", chr), .(g, og, chr)],
            focus_auto[, .(og, achr)], by = "og")
sa <- merge(pg[genome == SOMA & flag %chin% CONF & grepl("^OZ", chr), .(g, og, chr)],
            focus_auto[, .(og, achr)], by = "og")
jmap <- ja[, .(N = uniqueN(g)), by = .(achr, chr)][order(achr, -N)][!duplicated(achr),
            .(achr, j_home_chr = chr)]
smap <- sa[, .(N = uniqueN(g)), by = .(achr, chr)][order(achr, -N)][!duplicated(achr),
            .(achr, s_home_chr = chr)]
focus_auto <- merge(focus_auto, jmap, by = "achr", all.x = TRUE)
focus_auto <- merge(focus_auto, smap, by = "achr", all.x = TRUE)

# per og: juncea copies on the expected home (j_home) and Soma copies on the expected Soma home (s)
jhome <- merge(pg[genome == AJ & flag %chin% CONF, .(g, og, chr)],
               focus_auto[, .(og, j_home_chr)], by = "og")[
                 chr == j_home_chr, .(j_home = uniqueN(g)), by = og]
shome <- merge(pg[genome == SOMA & flag %chin% CONF, .(g, og, chr)],
               focus_auto[, .(og, s_home_chr)], by = "og")[
                 chr == s_home_chr, .(s = uniqueN(g)), by = og]

G_auto <- focus_auto
for (x in list(jhome, shome, jcnt, og_nchr)) G_auto <- merge(G_auto, x, by = "og", all.x = TRUE)
for (cc in c("j_home","s","neoX","ancX","Y","auto","unplaced","nchr"))
  G_auto[is.na(get(cc)), (cc) := 0L]
G_auto[, keep := nchr <= OVERMERGE_CHR & og %chin% sc]

cat(sprintf("\n== denominators (keep = nchr <= %d & strict 1:1:1) ==\n  neo-X: %d   anc-X: %d   auto: %d\n",
            OVERMERGE_CHR,
            G_neoX[keep == TRUE, .N], G_ancX[keep == TRUE, .N], G_auto[keep == TRUE, .N]))

wilson <- function(k, n, z = 1.96) {
  if (is.na(n) || n == 0) return(c(NA, NA))
  p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d; hw <- z * sqrt(p*(1-p)/n + z^2/(4*n^2))/d
  c(max(0, ctr - hw), min(1, ctr + hw))
}

# one panel per destination
panels <- list(
  auto = list(dest = "auto", main = "moved to placed autosomes",   ylab = "% moved to autosomes"),
  neoX = list(dest = "neoX", main = "moved onto the neo-X arm",    ylab = "% moved onto neo-X arm"),
  ancX = list(dest = "ancX", main = "moved onto the anc-X arm",    ylab = "% moved onto anc-X arm"),
  Y    = list(dest = "Y",    main = "moved onto the neo-Y",         ylab = "% moved onto neo-Y"))

# movers onto each destination: kept og, Soma keeps a copy on the homolog, juncea has none on the home arm
# (nor on the sister arm for neo-X ancestry) and >= 1 on the destination. Self bars (own arm; neo-Y for
# neo-X ancestry) are structurally zero. Then rate, Wilson CI, and pairwise Fisher tests of non-self bars
results <- list()
fisher_rows <- list()
for (nm in names(panels)) {
  dest <- panels[[nm]]$dest
  self_neoX <- dest %in% c("neoX", "Y")
  self_ancX <- dest == "ancX"
  k_neoX <- if (self_neoX) 0L else
    G_neoX[keep == TRUE & neoX == 0L & Y == 0L & s >= 1L & get(dest) >= 1L, .N]
  k_ancX <- if (self_ancX) 0L else
    G_ancX[keep == TRUE & ancX == 0L & s >= 1L & get(dest) >= 1L, .N]
  # autosomal ancestry: for dest == "auto" this means moved to a DIFFERENT autosome
  k_auto <- G_auto[keep == TRUE & j_home == 0L & s >= 1L & get(dest) >= 1L, .N]
  d <- data.frame(
    class = c("neo-X", "ancestral X", "autosomes"),
    k     = c(k_neoX,  k_ancX,        k_auto),
    n     = c(G_neoX[keep == TRUE, .N], G_ancX[keep == TRUE, .N], G_auto[keep == TRUE, .N]),
    self  = c(self_neoX, self_ancX, FALSE),
    stringsAsFactors = FALSE)
  d$rate <- ifelse(d$n == 0, NA_real_, d$k / d$n)
  ci <- t(mapply(wilson, d$k, d$n)); d$lo <- ci[, 1]; d$hi <- ci[, 2]
  results[[nm]] <- d
  ns <- d[!d$self, ]
  if (nrow(ns) >= 2L) {
    for (i in seq_len(nrow(ns) - 1L)) for (j in (i + 1L):nrow(ns)) {
      tb <- rbind(c(ns$k[i], ns$n[i] - ns$k[i]),
                  c(ns$k[j], ns$n[j] - ns$k[j]))
      ft <- suppressWarnings(fisher.test(tb))
      fisher_rows[[length(fisher_rows) + 1]] <- data.frame(
        panel  = nm,
        class1 = ns$class[i],
        class2 = ns$class[j],
        OR     = unname(ft$estimate),
        p      = ft$p.value,
        stringsAsFactors = FALSE)
    }
  }
}
# BH-FDR across all tests of all panels
fr <- do.call(rbind, fisher_rows)
fr$p.bh <- p.adjust(fr$p, method = "BH")

# print each panel with its tests
cat(sprintf("\n== Fisher's exact tests (BH-FDR across %d tests) ==\n", nrow(fr)))
for (nm in names(panels)) {
  cat(sprintf("\n-- panel: %s --\n", panels[[nm]]$main))
  print(results[[nm]])
  sub <- fr[fr$panel == nm, ]
  if (nrow(sub)) {
    for (r in seq_len(nrow(sub))) {
      cat(sprintf("  %-14s vs %-14s  OR=%.2f  p=%.4f  p.BH=%.4f\n",
                  sub$class1[r], sub$class2[r], sub$OR[r], sub$p[r], sub$p.bh[r]))
    }
  }
}

# anc-X -> neo-X movers: distance from the breakpoint (for the Results caveat)
mv <- G_ancX[keep == TRUE & ancX == 0L & s >= 1L & neoX >= 1L, og]
if (length(mv)) {
  mvd <- unique(pg[genome == AJ & og %chin% mv & chr == AJ_X, .(og, g, start)])
  mvd[, kb_from_brk := (start - BRK) / 1e3]
  cat("\n[anc-X -> neo-X movers: distance from breakpoint]\n")
  print(mvd[order(start)])
}

# fig4.R-style 'emigrated' (any off-pair destination, Y gametolog allowed) for comparison with "moved"
emig_neoX_fate <- G_neoX[keep == TRUE & neoX == 0L & (ancX + auto + unplaced) >= 1L & s >= 1L, .N]
emig_ancX_fate <- G_ancX[keep == TRUE & ancX == 0L & (neoX + auto + Y + unplaced) >= 1L & s >= 1L, .N]
cat(sprintf("\n[reconcile] fig4.R-style 'emigrated' (any off-pair destination):\n"))
cat(sprintf("  neo-X: %d / %d      anc-X: %d / %d\n",
            emig_neoX_fate, G_neoX[keep == TRUE, .N],
            emig_ancX_fate, G_ancX[keep == TRUE, .N]))

# horizontal bars as in fig4.R, bottom -> top: autosomes, neo-X, ancestral X. No device is opened,
# so the four panels go to the current device (Rplots.pdf under Rscript)
chrom_col <- c("neo-X" = "#7570B3", "ancestral X" = "#66A61E", "autosomes" = "#7F7F7F")
for (nm in names(panels)) {
  d <- results[[nm]][match(c("autosomes", "neo-X", "ancestral X"), results[[nm]]$class), ]
  d$rate <- ifelse(d$n == 0, 0, d$k / d$n)
  ci <- t(mapply(wilson, d$k, d$n)); d$lo <- ci[, 1]; d$hi <- ci[, 2]
  xmax <- min((max(c(d$hi, d$rate), na.rm = TRUE) * 1.35 + 0.005), 1)
  if (!is.finite(xmax) || xmax <= 0) xmax <- 0.05

  op <- par(mar = c(4.0, 5.5, 2.6, 4.0), las = 1, mgp = c(3.0, 0.7, 0))
  bp <- barplot(d$rate, names.arg = d$class, col = chrom_col[d$class], border = NA,
                xlim = c(0, xmax), horiz = TRUE, xlab = "proportion of ancestral genes moved",
                main = panels[[nm]]$main, cex.names = 1.0, cex.lab = 1.0, cex.main = 1.02)
  cap <- 0.15
  segments(d$lo, bp, d$hi, bp, lwd = 1.5, col = "grey25")
  segments(d$lo, bp - cap, d$lo, bp + cap, lwd = 1.5, col = "grey25")
  segments(d$hi, bp - cap, d$hi, bp + cap, lwd = 1.5, col = "grey25")
  text(pmax(d$hi, d$rate, na.rm = TRUE) + 0.02 * xmax, bp,
       sprintf("%.2g%% (%d/%d)", d$rate * 100, d$k, d$n),
       adj = c(0, 0.5), cex = 0.8, xpd = NA)
  if (any(d$self)) {
    mtext("(home / fusion homolog)", side = 1, line = 2.2, at = bp[d$self],
          cex = 0.72, col = "grey45")
  }
  abline(v = 0); par(op)
}
