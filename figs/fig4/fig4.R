# Fig. 4: copy-number fate of ancestral Anax genes on the A. juncea autosomes, anc-X, neo-X and neo-Y
suppressMessages(library(data.table))
PG   <- "../../odonate_synteny/genespace/pangenes/Aeshna_juncea_pangenes.txt.gz"
AJ   <- "Aeshna_juncea";        AJ_X <- "OZ369554.1"; AJ_Y <- "OZ369555.1"
ANAX <- "Anax_junius";          ANAX_NEOX <- "OZ349811.1"; ANAX_ANCX <- "OZ349813.1"
SOMA <- "Somatochlora_arctica"; SOMA_NEOX <- "OZ407514.1"; SOMA_ANCX <- "OZ407525.1"
CONF <- c("PASS", "array"); OVERMERGE_CHR <- 3
base <- "../../gene_pairs_ds"; P_ <- function(...) file.path(base, ...)   # pairs.sbatch run directory
CMP_GAMETO <- "juncea_X__vs__juncea_Y"   # reciprocal neoX<->neoY gametolog anchors

# anchor-partner rank filter shared with figs/fig3/fig3c.R; keep these thresholds equal to the ones there
source("../../lib/anchor_filter.R")
min.ratio     <- 0.65   # rank floor
min.pid       <- 0      # identity floor disabled: the rule is rank alone
ratio.na.keep <- FALSE

# fusion breakpoint written by gene_pairs_ds/pairs.sbatch (103 Mb if the file is missing)
BP_FILE <- P_("breakpoint.txt")
BRK <- if (file.exists(BP_FILE)) {
  as.numeric(readLines(BP_FILE))
} else {
  warning("breakpoint.txt not found; falling back to 103e6")
  103e6
}
cat(sprintf("[breakpoint] using %.3f Mb (from %s)\n", BRK / 1e6, BP_FILE))
gid  <- function(x) sub("\\.t[0-9]+$", "", x)

# TOGA2 join outputs (05_join.R), copied from gene_fate/toga2_out/join/ on the cluster
TOGA_DIR  <- "../../gene_fate/toga2_out/join"
DROP_FILE <- file.path(TOGA_DIR, "juncea_broken_genes.txt")   # set to NULL for the unfiltered run
USE_TOGA_RESCUE <- TRUE   # count intact TOGA2 copies missed by the annotation (tables print both ways)

# juncea genes with a TOGA2-verified broken CDS: dropped from copy counts and anchor pairs
drop <- if (!is.null(DROP_FILE) && file.exists(DROP_FILE)) readLines(DROP_FILE) else character()
if (!is.null(DROP_FILE) && !file.exists(DROP_FILE)) warning("[toga] ", DROP_FILE, " not found -- filter OFF")
cat(sprintf("[toga] %d juncea genes on the brokenCDS drop list%s\n", length(drop),
            if (!length(drop)) " (filter OFF)" else ""))

# reciprocal-singleton MCScanX anchors of one comparison, rank-filtered (transcript ids until the end)
synteny_recip_pairs <- function(cmp) {
  parts <- strsplit(cmp, "__vs__")[[1]]
  mp <- P_("comp", cmp, "mcscanx.pairs.tsv"); ida <- P_("chrom", paste0(parts[1], ".ids"))
  if (!file.exists(mp) || file.size(mp) == 0 || !file.exists(ida)) {
    warning(sprintf("[%s] no anchors", cmp)); return(data.table(gA = character(), gB = character()))
  }
  m <- fread(mp, header = FALSE)[, 1:2]; setnames(m, c("g1", "g2")); aids <- readLines(ida)
  n_pre <- nrow(m); m <- m[!(gid(g1) %chin% drop | gid(g2) %chin% drop)]
  cat(sprintf("[%s] %d anchor pairs, %d removed (brokenCDS partner)\n", cmp, n_pre, n_pre - nrow(m)))
  a1 <- m$g1 %in% aids
  d <- unique(data.table(gene_a = fifelse(a1, m$g1, m$g2),
                         gene_b = fifelse(a1, m$g2, m$g1)))
  d[, `:=`(kA = gid(gene_a), kB = gid(gene_b))]
  d <- d[, nA := .N, by = kA][, nB := .N, by = kB][nA == 1L & nB == 1L]
  if (!nrow(d)) return(data.table(gA = character(), gB = character()))

  n0 <- nrow(d)
  d  <- attach_rank(d, cmp, base)
  cat(sprintf("[%s] anchors %d -> %d after rank filter (min.ratio = %.2f)\n",
              cmp, n0, nrow(d), min.ratio))
  rank_report(d, cmp)
  unique(data.table(gA = gid(d$gene_a), gB = gid(d$gene_b)))
}

# pangenes; juncea genes on the drop list are re-flagged brokenCDS (non-confident)
pg <- fread(PG); pg[, g := gid(id)]; pg[, start := as.numeric(start)]; pg[is.na(og), og := -pgID]
n_flag <- pg[genome == AJ & g %chin% drop & flag %chin% CONF, .N]
pg[genome == AJ & g %chin% drop, flag := "brokenCDS"]
cat(sprintf("[toga] %d confident juncea pangenes rows re-flagged brokenCDS\n", n_flag))
if (length(drop) && n_flag == 0)
  warning("[toga] drop list loaded but no pangenes row matched: juncea ids in the pangenes file (e.g. ",
          pg[genome == AJ, id][1], ") do not look like the drop list (e.g. ", drop[1], ")")

# juncea arm of each gene, confident juncea copies per og, and Anax chromosomes per og (over-merge filter)
jloc <- unique(pg[genome == AJ, .(g, chr, start)])
jloc[, arm := fifelse(chr == AJ_X & start <  BRK, "neoX",
              fifelse(chr == AJ_X & start >= BRK, "ancX",
              fifelse(chr == AJ_Y,                "Y",
              fifelse(grepl("^OZ", chr),          "auto", "other"))))]
jloc <- unique(jloc[, .(jg = g, arm)])
jc      <- merge(pg[genome == AJ & flag %chin% CONF, .(g, og)], jloc, by.x = "g", by.y = "jg")
og_nchr <- pg[genome == ANAX & grepl("^OZ", chr), .(nchr = uniqueN(chr)), by = og]

# gametolog rescue: ogs whose neo-X member has a reciprocal neo-Y anchor (Y copy split into another og)
y_surv     <- unique(synteny_recip_pairs(CMP_GAMETO)$gA)       # neo-X genes with a surviving neo-Y gametolog
og_yanchor <- jc[arm == "neoX" & g %chin% y_surv, unique(og)]  # ogs anchor-present on the neo-Y

# TOGA2 presence: Anax genes projected intact (FI/I), orthologous and not retro / intronless, per juncea arm
TPROJ <- file.path(TOGA_DIR, "anax_projections.tsv")
tproj <- if (file.exists(TPROJ)) fread(TPROJ, colClasses = list(character = "juncea")) else NULL
if (is.null(tproj)) {
  if (USE_TOGA_RESCUE) warning("[toga] ", TPROJ, " not found (rerun the current 05_join.R) -- TOGA2 rescue OFF")
  USE_TOGA_RESCUE <- FALSE
  toga_on <- data.table(ag = character(), arm = character())
} else {
  tproj[is.na(juncea), juncea := ""]
  tproj[, retroish := retro | (exons == 1L & !is.na(ref_exons) & ref_exons >= 2L)]
  toga_on <- unique(tproj[status %chin% c("FI", "I") & !para & !retroish, .(ag, arm = region)])
  cat("[toga] intact orthologous non-retro TOGA2 copies (Anax genes) by juncea arm:\n")
  print(toga_on[, .N, by = arm][order(arm)])
  cat(sprintf("[toga] intact orthologous copies set aside as retro / intronless: %d projections\n",
              tproj[status %chin% c("FI", "I") & !para & retroish, .N]))
}

# Soma OZ407514.1 is neo-X homolog + auto-S, interleaved with no seam: a Soma gene counts as neo-X
# homolog only if its own og links to Anax neo-X (genes without a confident Anax link are missed)
og_achr <- pg[genome == ANAX & grepl("^OZ", chr), .N, by = .(og, chr)][order(-N), .(achr = chr[1]), by = og]
sf <- merge(pg[genome == SOMA & chr == SOMA_NEOX, .(g, og)], og_achr, by = "og", all.x = TRUE)
sf[, blk := fifelse(!is.na(achr) & achr == ANAX_NEOX, "neoXhom", "autoS")]
soma_neox_genes <- sf[blk == "neoXhom", g]                    # neo-X-homolog Soma genes (by Anax orthology)
soma_ancx_genes <- pg[genome == SOMA & chr == SOMA_ANCX, g]   # anc-X homolog is unfused -> all genes
cat(sprintf("[Soma fusion split] OZ407514.1 by per-gene Anax linkage: neo-X-homolog %d genes, auto-S excluded %d\n",
            length(soma_neox_genes), nrow(sf) - length(soma_neox_genes)))

# per focal Anax gene: a, s = Anax / Soma copies on the homolog, t = juncea total, h = home arm,
# sis = sister arm, o = off the proto-sex pair; optional anchor rescue and TOGA2 rescue of h and sis
build_cn <- function(anax_src, soma_genes, home_arm, sister_arm, yanchor_og = NULL, toga = FALSE) {
  foc  <- pg[genome == ANAX & chr == anax_src & flag %chin% CONF, .(og = og[1]), by = .(ag = g)]
  acnt <- pg[genome == ANAX & chr == anax_src & flag %chin% CONF, .(a = uniqueN(g)), by = og]
  scnt <- pg[genome == SOMA & g %chin% soma_genes & flag %chin% CONF, .(s = uniqueN(g)), by = og]
  jtot <- jc[,                .(t = uniqueN(g)), by = og]
  jhom <- jc[arm == home_arm, .(h = uniqueN(g)), by = og]
  G <- foc
  for (x in list(acnt, scnt, jtot, jhom)) G <- merge(G, x, by = "og", all.x = TRUE)
  if (!is.na(sister_arm))
    G <- merge(G, jc[arm == sister_arm, .(sis = uniqueN(g)), by = og], by = "og", all.x = TRUE)
  G <- merge(G, og_nchr, by = "og", all.x = TRUE)
  if (is.na(sister_arm)) G[, sis := 0L]
  for (cc in c("a", "s", "t", "h", "sis", "nchr")) G[is.na(get(cc)), (cc) := 0L]
  G[, o := pmax(0L, t - h - sis)]                       # off-pair copies (from the orthogroup)
  G[, rescued := FALSE]
  if (!is.null(yanchor_og)) {                           # anchor rescue (neo-Y bar only)
    G[(og %in% yanchor_og) & h == 0L, rescued := TRUE]
    G[rescued == TRUE, h := 1L]                         # recover an og-split Y gametolog as 1 home-arm copy
    G[, t := h + sis + o]                               # recompute total consistently
  }
  # TOGA2 rescue: h or sis of 0 becomes the number of the og's focal genes TOGA2 finds intact on that arm
  G[, `:=`(tres_h = FALSE, tres_sis = FALSE)]
  if (toga && nrow(toga_on)) {
    for (side in c("h", "sis")) {
      a_ <- if (side == "h") home_arm else sister_arm
      if (is.na(a_)) next
      hit <- toga_on$ag[toga_on$arm == a_]
      th  <- G[ag %chin% hit, .(th = uniqueN(ag)), by = og]
      G   <- merge(G, th, by = "og", all.x = TRUE)
      G[is.na(th), th := 0L]
      G[get(side) == 0L & th > 0L, c(side, paste0("tres_", side)) := .(th, TRUE)]
      G[, th := NULL]
    }
    G[, t := h + sis + o]
  }
  G[]
}

# categories, first match wins; excluded = Soma-silent, over-merged og, or a != s (ancestral CN ambiguous)
classify_cn <- function(G) {
  G[, cat := fcase(
    t == 0 & s >= 1 & a >= 1,                           "fully deleted",
    h == 0 & o >= 1 & s >= 1,                           "emigrated",
    h == 0 & sis >= 1 & o == 0 & s >= 1,                "degenerated",
    h > pmax(a, s) & o == 0,                            "duplicated",
    pmin(a, s) >= 2 & h > 0 & h < pmin(a, s) & o == 0,  "partial deletion",
    h == a & a == s & o == 0,                           "conserved",
    a >= 1 & s >= 1,                                    "other",
    default = "excluded")]
  G[s == 0 | nchr > OVERMERGE_CHR | a != s, cat := "excluded"][]
}

# classify neo-X, anc-X and neo-Y, with and without the TOGA2 rescue
G_neoX <- classify_cn(build_cn(ANAX_NEOX, soma_neox_genes, "neoX", "Y",  toga = USE_TOGA_RESCUE))
G_ancX <- classify_cn(build_cn(ANAX_ANCX, soma_ancx_genes, "ancX", NA,  toga = USE_TOGA_RESCUE))
G_Y    <- classify_cn(build_cn(ANAX_NEOX, soma_neox_genes, "Y",    "neoX", yanchor_og = og_yanchor,
                               toga = USE_TOGA_RESCUE))
G0_neoX <- classify_cn(build_cn(ANAX_NEOX, soma_neox_genes, "neoX", "Y"))
G0_ancX <- classify_cn(build_cn(ANAX_ANCX, soma_ancx_genes, "ancX", NA))
G0_Y    <- classify_cn(build_cn(ANAX_NEOX, soma_neox_genes, "Y",    "neoX", yanchor_og = og_yanchor))

# pooled autosome baseline: home in juncea and Soma is the og's modal autosome (ties: first in sort order)
build_cn_auto <- function(toga = FALSE) {
  anax_auto <- setdiff(unique(pg[genome == ANAX & grepl("^OZ", chr), chr]),
                       c(ANAX_NEOX, ANAX_ANCX))
  soma_auto <- setdiff(unique(pg[genome == SOMA & grepl("^OZ", chr), chr]),
                       c(SOMA_NEOX, SOMA_ANCX))
  # focal: one row per Anax gene on any Anax autosome; a = confident copies on the same Anax chromosome
  foc <- pg[genome == ANAX & chr %chin% anax_auto & flag %chin% CONF,
            .(og = og[1], anax_src = chr[1]), by = .(ag = g)]
  acnt <- pg[genome == ANAX & chr %chin% anax_auto & flag %chin% CONF,
             .(a = uniqueN(g)), by = .(og, anax_src = chr)]
  # h = juncea copies on the og's modal juncea autosome
  jc_auto <- pg[genome == AJ & flag %chin% CONF & grepl("^OZ", chr) &
                !chr %chin% c(AJ_X, AJ_Y), .(g, og, chr)]
  jc_mode <- jc_auto[, .N, by = .(og, chr)][order(og, -N)][
             , .SD[1], by = og][, .(og, jhome_chr = chr)]
  jhom <- merge(jc_auto, jc_mode, by = "og")[chr == jhome_chr,
                                              .(h = uniqueN(g)), by = og]
  # s = Soma copies on the og's modal Soma autosome
  soma_at   <- pg[genome == SOMA & flag %chin% CONF & chr %chin% soma_auto,
                  .(g, og, chr)]
  soma_mode <- soma_at[, .N, by = .(og, chr)][order(og, -N)][
               , .SD[1], by = og][, .(og, shome_chr = chr)]
  scnt <- merge(soma_at, soma_mode, by = "og")[chr == shome_chr,
                                                .(s = uniqueN(g)), by = og]
  # t = total confident juncea copies genome-wide
  jtot <- jc[, .(t = uniqueN(g)), by = og]

  G <- merge(foc, acnt, by = c("og", "anax_src"), all.x = TRUE)
  for (x in list(jtot, jhom, scnt, og_nchr)) G <- merge(G, x, by = "og", all.x = TRUE)
  G[, sis := 0L]                            # no sister-arm concept for autosomes
  for (cc in c("a", "s", "t", "h", "sis", "nchr"))
    G[is.na(get(cc)), (cc) := 0L]
  G[, o := pmax(0L, t - h - sis)]
  G[, rescued := FALSE]                     # no anchor rescue for autosomes
  # TOGA2 rescue: if t == 0, h = focal genes TOGA2 finds intact on any juncea autosome
  G[, `:=`(tres_h = FALSE, tres_sis = FALSE)]
  if (toga && nrow(toga_on)) {
    hit <- toga_on$ag[toga_on$arm == "auto"]
    th  <- G[ag %chin% hit, .(th = uniqueN(ag)), by = og]
    G   <- merge(G, th, by = "og", all.x = TRUE)
    G[is.na(th), th := 0L]
    G[t == 0L & th > 0L, `:=`(h = th, tres_h = TRUE)]
    G[, th := NULL][, t := h + sis + o]
  }
  G[]
}
G_auto  <- classify_cn(build_cn_auto(toga = USE_TOGA_RESCUE))
G0_auto <- classify_cn(build_cn_auto())

# the four chromosome classes under each classification
arms   <- c("autosomes", "ancestral X", "neo-X", "neo-Y")
G_sets <- list(toga2 = setNames(list(G_auto,  G_ancX,  G_neoX,  G_Y),  arms),
               annot = setNames(list(G0_auto, G0_ancX, G0_neoX, G0_Y), arms))

cat("\n== pangenes rows vs distinct genes (if these differ, .N and uniqueN(g) disagree) ==\n")
print(pg[, .(rows = .N, genes = uniqueN(g)), by = genome])

# per class: what each exclusion removes (anc_mismatch = cost of requiring a == s), what is left in "other"
dr <- list()
for (nm in arms) {
  G <- G_sets$toga2[[nm]]
  dr[[nm]] <- data.table(arm          = nm,
                         focal        = nrow(G),
                         soma_silent  = G[s == 0, .N],
                         overmerged   = G[s > 0 & nchr > OVERMERGE_CHR, .N],
                         anc_mismatch = G[s > 0 & nchr <= OVERMERGE_CHR & a != s, .N],
                         informative  = G[!(s == 0 | nchr > OVERMERGE_CHR | a != s), .N],
                         other_left   = G[cat == "other", .N])
}
cat("\n== exclusions, and the residual that survives them ==\n")
print(rbindlist(dr))
cat("  anc_mismatch is the cost of requiring a == s. Compare its rate across\n",
    "  arms: a much higher rate on autosomes than on the X arms would bias the\n",
    "  autosomal baseline upward and widen every sex-chromosome contrast.\n",
    "  other_left should now be only the h >= 1 & o >= 1 case.\n", sep = "")

# category counts per class, in copy-number order, then the non-CN categories and the excluded set
segs <- c("fully deleted", "partial deletion", "conserved", "duplicated",
          "degenerated", "emigrated", "other", "excluded")
tabf <- function(G) as.integer(table(factor(G$cat, segs)))
M <- rbind(`autosomes`   = tabf(G_auto),
           `ancestral X` = tabf(G_ancX),
           `neo-X`       = tabf(G_neoX),
           `neo-Y`       = tabf(G_Y))
colnames(M) <- segs

# the excluded set must be the same for neo-X and neo-Y, and the categories must partition each focal set
stopifnot(M["neo-X", "excluded"] == M["neo-Y", "excluded"])
stopifnot(all(rowSums(M) == c(nrow(G_auto), nrow(G_ancX), nrow(G_neoX), nrow(G_Y))))

# denominator = informative set; then two composite "reduced" categories (with and without full loss)
M_all <- M                                                # keep the full table incl. excluded
M     <- M[, setdiff(segs, "excluded"), drop = FALSE]
tot   <- rowSums(M)                                       # denominator computed BEFORE composites
M <- cbind(M,
           `reduced (incl. fully lost)` = M[, "partial deletion"] + M[, "fully deleted"],
           `reduced (excl. fully lost)` = M[, "partial deletion"])

cat(sprintf("\nfocal Anax genes  |  autosomes: %d   anc-X: %d   neo-X: %d   neo-Y: %d\n",
            nrow(G_auto), nrow(G_ancX), nrow(G_neoX), nrow(G_Y)))
cat(sprintf("excluded          |  autosomes: %d   anc-X: %d   neo-X: %d   neo-Y: %d\n",
            M_all["autosomes",   "excluded"], M_all["ancestral X", "excluded"],
            M_all["neo-X",       "excluded"], M_all["neo-Y",       "excluded"]))
cat("\n== counts (informative set: excluded ogs dropped) ==\n"); print(M)
cat("\n== proportions ==\n"); print(round(M / tot, 3))
cat(sprintf("\n[gametolog rescue] neo-Y genes recovered by anchor set: %d\n",
            G_Y[rescued == TRUE, .N]))
cat(sprintf("[TOGA2 rescue %s] genes given a home-arm copy: autosomes %d  anc-X %d  neo-X %d  neo-Y %d;  sister-arm copy: neo-X %d  neo-Y %d\n",
            if (USE_TOGA_RESCUE) "ON" else "OFF",
            G_auto[tres_h == TRUE, .N], G_ancX[tres_h == TRUE, .N], G_neoX[tres_h == TRUE, .N],
            G_Y[tres_h == TRUE, .N], G_neoX[tres_sis == TRUE, .N], G_Y[tres_sis == TRUE, .N]))

# category counts without -> with the TOGA2 rescue
cmp <- NULL
for (nm in arms) {
  n0 <- tabf(G_sets$annot[[nm]]); n1 <- tabf(G_sets$toga2[[nm]])
  cmp <- rbind(cmp, ifelse(n0 == n1, as.character(n0), paste0(n0, "->", n1)))
}
dimnames(cmp) <- list(arms, segs)
cat("\n== counts without -> with the TOGA2 rescue (all categories incl. excluded) ==\n")
print(noquote(cmp))

# TOGA2-only presence (sensitivity): per informative focal gene, is it projected intact on its arm,
# only on the sister arm, only elsewhere, or nowhere; presence only, no copy number
if (nrow(toga_on)) {
  ARMS <- c("Y", "neoX", "ancX", "auto", "other")
  tog_only <- function(G, home, sister = NA) {
    x   <- G[cat != "excluded", .(ag, cat)]
    on  <- function(arms) x$ag %chin% toga_on[arm %chin% arms, ag]
    x[, `:=`(H = on(home), S = if (is.na(sister)) FALSE else on(sister),
             E = on(setdiff(ARMS, c(home, sister))))]
    x[, tcat := fcase(H, "intact on arm", S, "intact on sister only",
                      E, "intact elsewhere only", default = "no intact copy")][]
  }
  T_auto <- tog_only(G_auto, "auto"); T_ancX <- tog_only(G_ancX, "ancX")
  T_neoX <- tog_only(G_neoX, "neoX", "Y"); T_Y <- tog_only(G_Y, "Y", "neoX")
  tl <- c("intact on arm", "intact on sister only", "intact elsewhere only", "no intact copy")
  tt <- rbind(autosomes = table(factor(T_auto$tcat, tl)), `ancestral X` = table(factor(T_ancX$tcat, tl)),
              `neo-X` = table(factor(T_neoX$tcat, tl)), `neo-Y` = table(factor(T_Y$tcat, tl)))
  cat("\n== TOGA2-only presence (informative set; 'intact on sister only' = degenerated) ==\n")
  print(cbind(tt, n = rowSums(tt)))
  print(round(tt / rowSums(tt), 3))

  # neo arms side by side: annotation only, TOGA2 rescue, TOGA2 only
  sbs <- function(lab, G0, G1, Tn) {
    pc <- function(v, n) sprintf("%5.1f%%", 100 * v / n)
    n <- nrow(Tn)
    data.table(arm = lab, n = n,
      `conserved: annot` = pc(G0[cat == "conserved", .N], n), `rescue` = pc(G1[cat == "conserved", .N], n),
      `toga2-only (intact on arm)` = pc(Tn[tcat == "intact on arm", .N], n),
      `degenerated: annot` = pc(G0[cat == "degenerated", .N], n), `rescue ` = pc(G1[cat == "degenerated", .N], n),
      `toga2-only (sister only)` = pc(Tn[tcat == "intact on sister only", .N], n))
  }
  cat("\n== neo arms: annotation only / TOGA2 rescue / TOGA2 only ==\n")
  print(rbind(sbs("neo-X", G0_neoX, G_neoX, T_neoX), sbs("neo-Y", G0_Y, G_Y, T_Y)))
  cat("   (conserved in the main table also requires h == a; TOGA2-only is presence, so compare loosely)\n")
  cat("\n== neo-Y: rescued classification x TOGA2-only ==\n")
  print(dcast(T_Y[, .N, by = .(cat, tcat = factor(tcat, tl))], cat ~ tcat, value.var = "N", fill = 0L, drop = FALSE))
  fwrite(rbindlist(list(T_auto[, arm := "autosomes"], T_ancX[, arm := "ancestral X"],
                        T_neoX[, arm := "neo-X"], T_Y[, arm := "neo-Y"])), "toga2_only_calls.tsv", sep = "\t")
}

wilson <- function(k, n, z = 1.96) {
  if (is.na(n) || n == 0) return(c(NA, NA))
  p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d; hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d
  c(max(0, ctr-hw), min(1, ctr+hw))
}

# proportion of each class's ancestral genes per category, with Wilson 95% CI
vals <- list()
for (cc in setdiff(segs, "excluded")) {
  k <- M[, cc]; ci <- t(mapply(wilson, k, tot))
  vals[[cc]] <- data.table(category = cc, chrom = rownames(M), k = as.integer(k), n = as.integer(tot),
                           prop = round(k/tot, 4), lo = round(ci[, 1], 4), hi = round(ci[, 2], 4))
}
vals <- rbindlist(vals)
cat("\n== values (proportion of each chromosome's ancestral genes, Wilson 95% CI) ==\n"); print(vals)

# figure colours and titles
chrom_col <- c("neo-X"       = "#7570B3",
               "ancestral X" = "#66A61E",
               "neo-Y"       = "#B3717C",
               "autosomes"   = "#888888")
ctitle <- c("fully deleted"     = "lost from genome",
            "partial deletion"  = "reduced copy number on chromosome",
            "conserved"         = "copy number conserved on chromosome",
            "duplicated"        = "amplified copy number on chromosome",
            "degenerated"       = "lost from chromosome, retained on sister arm",
            "emigrated"         = "relocated off the sex-chromosome pair",
            "other"             = "ancestrally present, non-CN-spectrum classification",
            "reduced (incl. fully lost)" = "reduced copy number (partial + full loss)",
            "reduced (excl. fully lost)" = "reduced copy number (partial only)")

# plotting matrices (built as M / tot above) for both classifications; identical if USE_TOGA_RESCUE is FALSE
set_label <- c(toga2 = "with TOGA2 rescue", annot = "annotation only")
fig_sets <- list()
for (tag in names(G_sets)) {
  m <- NULL
  for (G in G_sets[[tag]]) m <- rbind(m, tabf(G))
  dimnames(m) <- list(arms, segs)
  m <- m[, setdiff(segs, "excluded"), drop = FALSE]
  n <- rowSums(m)
  m <- cbind(m, `reduced (incl. fully lost)` = m[, "partial deletion"] + m[, "fully deleted"],
                `reduced (excl. fully lost)` = m[, "partial deletion"])
  fig_sets[[tag]] <- list(M = m, tot = n, label = set_label[[tag]])
}
stopifnot(identical(fig_sets$toga2$M, M), identical(fig_sets$toga2$tot, tot))
ord_chr <- c("neo-Y", "neo-X", "ancestral X", "autosomes")   # bars bottom -> top, autosomes on top

# 10 figures: 5 categories x 2 classifications; "degenerated" is 0 by construction on autosomes and anc-X
plot_cats <- c("conserved", "duplicated",
               "reduced (incl. fully lost)", "reduced (excl. fully lost)",
               "degenerated")
graphics.off()                                   # clear any stale devices first
outdir <- "."                                    # e.g. "figs" to keep them separate
dir.create(outdir, showWarnings = FALSE)

for (cc in plot_cats) {
  # shared x-axis for the two versions: largest upper CI or proportion, +30%
  xm <- c()
  for (S in fig_sets) {
    k <- S$M[ord_chr, cc]; n <- S$tot[ord_chr]; ci <- t(mapply(wilson, k, n))
    xm <- c(xm, min(max(c(ci[, 2], k / n), na.rm = TRUE) * 1.30 + 0.02, 1))
  }
  xm <- max(xm)

  for (tag in names(fig_sets)) {
    S  <- fig_sets[[tag]]
    f  <- file.path(outdir, paste0("fate.fig.", gsub("[^A-Za-z0-9]+", "_", gsub("[()]", "", cc)),
                                   ".", tag, ".svg"))
    # horizontal bars with Wilson CIs and "% (k/n)" labels
    ok <- tryCatch({
      svg(f, width = 18, height = 10)
      k <- S$M[ord_chr, cc]
      n <- S$tot[ord_chr]
      prop <- k / n
      ci <- t(mapply(wilson, k, n)); lo <- ci[, 1]; hi <- ci[, 2]
      op <- par(mar = c(4.0, 5.5, 2.6, 4.0), las = 1, mgp = c(3.0, 0.7, 0))
      bp <- barplot(prop, names.arg = ord_chr, col = chrom_col[ord_chr], border = NA, xlim = c(0, xm),
                    horiz = TRUE, xlab = "proportion of the chromosome's ancestral genes",
                    main = paste0(ctitle[[cc]], " (", S$label, ")"),
                    tcl = -2, mgp = c(3, 3.5, 0), cex.names = 1.0, cex.lab = 1.0, cex.main = 1.02)
      cap <- 0.15
      segments(lo, bp, hi, bp, lwd = 1.5, col = "grey25")
      segments(lo, bp - cap, lo, bp + cap, lwd = 1.5, col = "grey25")
      segments(hi, bp - cap, hi, bp + cap, lwd = 1.5, col = "grey25")
      text(pmax(hi, prop, na.rm = TRUE) + 0.02 * xm, bp, sprintf("%.1f%% (%d/%d)", prop*100, k, n),
           adj = c(0, 0.5), cex = 0.8, xpd = NA)
      abline(v = 0); par(op)
      dev.off(); TRUE
    }, error = function(e) { message("[FAIL] ", cc, " / ", tag, ": ", conditionMessage(e)); FALSE })
    if (!ok && dev.cur() != 1L) dev.off()        # close the orphaned device
    cat(sprintf("%-28s %-6s %-6s %s\n", cc, tag, ok,
                if (file.exists(f)) paste(file.size(f), "bytes") else "NO FILE"))
  }
}
cat("\nwd:", getwd(), "\n"); print(list.files(outdir, pattern = "^fate\\.fig\\..*\\.(toga2|annot)\\.svg$"))

# TOGA2 check of the annotation-only neo-Y "degenerated" genes (reviewer request): TOGA2's neo-Y call
TG <- file.path(TOGA_DIR, "anax_gene_status.tsv")
if (!file.exists(TG)) message("[toga] ", TG, " not found; skipping the degenerated x TOGA2 table") else {
  tg  <- fread(TG)
  # 'inactivated' split by TOGA2's own call: L = lost, UL = uncertain loss
  fine <- function(call, best) fifelse(call == "inactivated", paste0("inactivated (", best, ")"), call)
  deg_chk <- function(G) merge(G[cat == "degenerated", .(ag, og)],
                               tg[, .(ag = gene, Y_call, X_call, Y_best)], by = "ag", all.x = TRUE)
  chk <- deg_chk(G0_Y)
  cat(sprintf("\n== the %d neo-Y 'degenerated' genes (annotation only) x TOGA2 neo-Y call ==\n", nrow(chk)))
  print(chk[, .N, by = .(Y_call = fine(Y_call, Y_best))][order(-N)])
  cat("   (X_call of the same genes, should be mostly intact:)\n")
  print(chk[, .N, by = X_call][order(-N)])

  # degenerated genes with an intact neo-Y projection: why no Y copy was counted. Each intact projection
  # gets the first class in cls_order; the first five are rescued by the TOGA2 rescue, the last two are not
  if (is.null(tproj)) message("[toga] no anax_projections.tsv; skipping the breakdown") else {
    tp <- tproj[region == "Y" & status %chin% c("FI", "I")]
    pj <- pg[genome == AJ, .(conf = any(flag %chin% CONF), cog = list(unique(og[flag %chin% CONF]))), by = g]
    cls_order <- c("same og (unexpected)", "different og", "not PASS/array", "not in pangenes",
                   "no gene model", "retrocopy", "paralog only")
    first_cls <- function(x) x[which.min(match(x, cls_order))]
    cls1 <- function(para, retroish, juncea, focal_og) {
      if (isTRUE(para)) return("paralog only")
      if (isTRUE(retroish)) return("retrocopy")
      if (!nzchar(juncea)) return("no gene model")
      # class of each juncea gene the projection overlaps
      cl <- character()
      for (jj in strsplit(juncea, ",", fixed = TRUE)[[1]]) {
        x <- pj[g == jj]
        cl <- c(cl, if (!nrow(x)) "not in pangenes"
                    else if (!x$conf) "not PASS/array"
                    else if (focal_og %in% x$cog[[1]]) "same og (unexpected)"
                    else "different og")
      }
      first_cls(cl)
    }
    n_int <- chk[Y_call == "intact", .N]
    di <- merge(chk[Y_call == "intact", .(ag, og)], tp, by = "ag", allow.cartesian = TRUE)
    di[, cls := mapply(cls1, para, retroish, juncea, og)]
    dg <- di[, .(cls = first_cls(cls)), by = .(ag, og)]
    cat(sprintf("\n== the %d with an intact neo-Y projection: why the annotation-based count saw no Y copy ==\n", n_int))
    print(dg[, .N, by = cls][order(match(cls, cls_order))])
    if (nrow(dg) != n_int)
      warning("[toga] ", n_int - nrow(dg), " intact-Y genes had no FI/I row in ", TPROJ)
    fwrite(merge(dg, di[, .(ag, proj_cls = cls, status, para, retro, exons, ref_exons, juncea, chr, start, end, base)],
                 by = "ag")[order(cls, ag)], "degenerated_Y_intact.tsv", sep = "\t")
    cat("   per-projection detail: degenerated_Y_intact.tsv\n")
  }
  if (USE_TOGA_RESCUE) {
    chk1 <- deg_chk(G_Y)
    cat(sprintf("\n== after the TOGA2 rescue: %d neo-Y 'degenerated' genes x TOGA2 neo-Y call ==\n", nrow(chk1)))
    print(chk1[, .N, by = .(Y_call = fine(Y_call, Y_best))][order(-N)])
  }

  # genes counted present (conserved / duplicated) x TOGA2's call for the same arm; nothing is reclassified.
  # Doubted copies at a similar rate on neo-X (the control) and neo-Y are alignment noise, not degeneration
  lv <- c("intact", "partially_intact", "inactivated (UL)", "inactivated (L)", "missing", "other",
          "no projection", "not in TOGA2")
  pres_sets <- list(`neo-Y annot`  = list(G0_Y,    "Y"), `neo-Y rescue` = list(G_Y,    "Y"),
                    `neo-X annot`  = list(G0_neoX, "X"), `neo-X rescue` = list(G_neoX, "X"))
  pres <- list()
  for (nm in names(pres_sets)) {
    G <- pres_sets[[nm]][[1]]; side <- pres_sets[[nm]][[2]]
    t2 <- tg[, c("gene", paste0(side, c("_call", "_best"))), with = FALSE]
    setnames(t2, c("ag", "call", "best"))
    x  <- merge(G[cat %chin% c("conserved", "duplicated"), .(ag, og)], t2, by = "ag", all.x = TRUE)
    x[, call := fine(call, best)][is.na(call), call := "not in TOGA2"][
      call %like% "^no_[XY]_copy$", call := "no projection"]
    pres[[nm]] <- x[, .N, by = call][, set := nm]
  }
  pres <- rbindlist(pres)
  pres <- dcast(pres, factor(call, lv) ~ factor(set, unique(set)), value.var = "N", fill = 0L)
  setnames(pres, 1, "TOGA2 call")
  cat("\n== genes counted present (conserved + duplicated) x TOGA2 call for the same arm ==\n")
  print(pres)
}

# query sets for the miniprot check (06_miniprot.slurm), one row per set x focal Anax gene x juncea query:
# degenerated / conserved_Y = annotation-only neo-Y calls, query = the neo-X copy;
# ancX_control = anc-X conserved, query = the anc-X copy (no neo-Y gametolog, so a good neo-Y hit is
# a false positive)
mp_sets <- list(degenerated  = list(G0_Y[cat == "degenerated"],  "neoX"),
                conserved_Y  = list(G0_Y[cat == "conserved"],    "neoX"),
                ancX_control = list(G0_ancX[cat == "conserved"], "ancX"))
ms <- list()
for (nm in names(mp_sets)) {
  G <- mp_sets[[nm]][[1]]; arm_q <- mp_sets[[nm]][[2]]
  x <- merge(G[, .(ag, og, cat_annot = cat)], jc[arm == arm_q, .(og, jg = g)],
             by = "og", allow.cartesian = TRUE)
  ms[[nm]] <- x[, set := nm][]
}
ms <- rbindlist(ms, use.names = TRUE)
ms <- merge(ms, G_Y[, .(ag, cat_rescue = cat)], by = "ag", all.x = TRUE)       # NA for anc-X genes
if (file.exists(TG)) ms <- merge(ms, fread(TG)[, .(ag = gene, Y_call, Y_best)], by = "ag", all.x = TRUE)
setcolorder(ms, c("set", "ag", "og", "jg"))
fwrite(ms[order(set, ag, jg)], "miniprot_sets.tsv", sep = "\t")
cat("\n== miniprot query sets (write: miniprot_sets.tsv) ==\n")
print(ms[, .(anax_genes = uniqueN(ag), juncea_queries = uniqueN(jg)), by = set])

# neo-X classes for the read-depth check (gene_fate/depth_check/02_build_table.R): ogs with one
# consistent neo-Y category and a single neo-X gene
yca   <- unique(G_Y[cat %chin% c("conserved", "degenerated"), .(og, ag, cat)])
agree <- yca[, .(cat = cat[1], ok = uniqueN(cat) == 1L), by = og]
yc    <- agree[ok == TRUE]
yc[, class := fifelse(cat == "conserved", "neoX_conserved", "neoX_lost")]
nx <- unique(jc[arm == "neoX", .(gene = g, og)])
nx[, n_nx := .N, by = og]
out <- merge(yc[, .(og, class)], nx[n_nx == 1L, .(og, gene)], by = "og")
fwrite(out[, .(gene, class)], "gene_classes.tsv", sep = "\t", col.names = FALSE)
writeLines(drop, "toga2_dropped_genes.txt")
cat(sprintf("\n[gene_classes] %d genes (conserved %d, lost %d)\n",
            nrow(out), out[class == "neoX_conserved", .N], out[class == "neoX_lost", .N]))
cat(sprintf("[gene_classes] orthogroups removed: %d conflicting category, %d with >1 neo-X gene, %d with no neo-X gene model\n",
            agree[ok == FALSE, .N], yc[og %in% nx[n_nx > 1L, og], .N], yc[!og %in% nx$og, .N]))
