# Scores the 06_miniprot.slurm alignments against the figs/fig4/fig4.R sets (called by 06_miniprot.slurm)
# usage: Rscript 06_score_miniprot.R <run dir> <A. juncea gene GFF>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); RUN <- a[1]; QGFF <- a[2]
MP <- file.path(RUN, "miniprot")
ms <- fread(file.path(MP, "miniprot_sets.tsv"))
qok <- sub("^>", "", grep("^>", readLines(file.path(MP, "query.faa")), value = TRUE))

# miniprot GFF: each alignment is a ##PAF line followed by its mRNA line
L   <- readLines(file.path(MP, "neoY.miniprot.gff"))
paf <- L[startsWith(L, "##PAF\t")]
mr  <- L[grepl("\tmRNA\t", L, fixed = TRUE)]
stopifnot(length(paf) == length(mr))
f   <- tstrsplit(sub("^##PAF\t", "", paf), "\t", fixed = TRUE, keep = 1:9)
tagi <- function(x, k) { v <- sub(paste0(".*\t", k, ":i:(-?[0-9]+).*"), "\\1", x); as.integer(ifelse(v == x, "0", v)) }
attr1 <- function(x, k) { v <- sub(paste0(".*[\t;]", k, "=([^;\t]+).*"), "\\1", x); v[v == x] <- NA; v }
al <- data.table(jg = f[[1]], qlen = as.integer(f[[2]]), qs = as.integer(f[[3]]), qe = as.integer(f[[4]]),
                 chr = f[[6]], start = as.integer(f[[8]]) + 1L, end = as.integer(f[[9]]),
                 fs = tagi(paf, "fs"), st = tagi(paf, "st"),
                 rank = as.integer(attr1(mr, "Rank")), ident = as.numeric(attr1(mr, "Identity")))
stopifnot(identical(sub(" .*", "", attr1(mr, "Target")), al$jg))   # PAF and mRNA lines pair up

# class of each query's best (Rank 1) alignment: intact >= 90% aligned, no frameshift or stop;
# broken >= 50% aligned with a frameshift or stop; partial otherwise
al[, cov := (qe - qs) / qlen]
al[, cls := fcase(fs + st > 0 & cov >= 0.5, "broken", fs + st == 0 & cov >= 0.9, "intact", default = "partial")]
al <- al[rank == 1L]

# TOGA2 projections on the neo-Y and A. juncea gene models, to check overlaps
al[, chr0 := sub("\\.[0-9]+$", "", chr)]
tp <- fread(file.path(RUN, "toga2_out/join/anax_projections.tsv"), colClasses = list(character = "juncea"))
tp <- tp[region == "Y", .(ag, tchr = chr, tstart = start, tend = end, tstatus = status, tpara = para)]
gf <- fread(cmd = paste("awk -F'\\t' '$3==\"gene\"' ", shQuote(QGFF)), header = FALSE, select = c(1, 4, 5, 9))
setnames(gf, c("chr", "gstart", "gend", "gid")); gf[, gid := sub(".*ID=([^;]+).*", "\\1", gid)]

# one row per set x Anax gene x query, then the best query per gene (intact > broken > partial > none)
lv <- c("intact", "broken", "partial", "none", "no query")
q <- merge(ms, al[, .(jg, chr, chr0, start, end, cov, fs, st, ident, cls)], by = "jg", all.x = TRUE)
q[is.na(cls), cls := fifelse(jg %chin% qok, "none", "no query")]
q[, `:=`(ord = match(cls, lv), ident0 = fifelse(is.na(ident), -1, ident))]
setorder(q, set, ag, ord, -ident0)
g <- q[, .SD[1], by = .(set, ag)]

# gametolog identity: 5th percentile of identity among intact positive-control hits
thr <- g[set == "conserved_Y" & cls == "intact", quantile(ident, 0.05, names = FALSE)]
g[, mp := factor(fifelse(cls == "intact" & ident >= thr, "intact, gametolog identity",
                 fifelse(cls == "intact", "intact, lower identity", cls)),
                 c("intact, gametolog identity", "intact, lower identity", lv[-1]))]

# same locus as a TOGA2 projection of the same Anax gene (any status), and overlap with a neo-Y gene model
ov <- merge(g[!is.na(start), .(set, ag, chr0, start, end)], tp, by = "ag", allow.cartesian = TRUE)[
  tchr == chr0 & tstart <= end & tend >= start,
  .(toga_locus = paste(sort(unique(tstatus)), collapse = ",")), by = .(set, ag)]
g <- merge(g, ov, by = c("set", "ag"), all.x = TRUE)
gm <- merge(g[!is.na(start), .(set, ag, chr, start, end)], gf, by = "chr", allow.cartesian = TRUE)[
  gstart <= end & gend >= start, .(y_model = paste(unique(gid), collapse = ",")), by = .(set, ag)]
g <- merge(g, gm, by = c("set", "ag"), all.x = TRUE)

# report
pc <- function(k, n) sprintf("%d (%.1f%%)", k, 100 * k / n)
cat(sprintf("\n[miniprot] gametolog identity threshold = %.3f (5th percentile of intact positive-control hits, n = %d)\n",
            thr, g[set == "conserved_Y" & cls == "intact", .N]))
cat("\n== best neo-Y alignment per gene, by set ==\n")
tab <- dcast(g[, .N, by = .(set, mp)], mp ~ set, value.var = "N", fill = 0L, drop = FALSE)
print(tab[, c("mp", intersect(c("conserved_Y", "degenerated", "ancX_control"), names(tab))), with = FALSE])
cat("\n== amino-acid identity of intact hits (quantiles) ==\n")
print(g[cls == "intact", as.list(round(quantile(ident, c(0, .05, .25, .5, .75, 1)), 3)), by = set])
fp <- g[set == "ancX_control", mean(mp == "intact, gametolog identity")]
cat(sprintf("false-positive rate (anc-X proteins intact at gametolog identity on the neo-Y): %.1f%%\n", 100 * fp))

# degenerated genes against the TOGA2 neo-Y call and the classification after the TOGA2 rescue
d <- g[set == "degenerated"]
if ("Y_call" %in% names(d)) {
  d[, toga := fifelse(Y_call == "inactivated", paste0("inactivated (", Y_best, ")"), Y_call)]
  cat("\n== degenerated (annotation only): miniprot x TOGA2 neo-Y call ==\n")
  print(dcast(d[, .N, by = .(mp, toga)], mp ~ toga, value.var = "N", fill = 0L, drop = c(TRUE, TRUE)))
}
if ("cat_rescue" %in% names(d)) {
  cat("\n== degenerated (annotation only): miniprot x classification after the TOGA2 rescue ==\n")
  print(dcast(d[, .N, by = .(mp, cat_rescue)], mp ~ cat_rescue, value.var = "N", fill = 0L))
  r <- d[cat_rescue == "conserved"]; k <- r[mp == "intact, gametolog identity", .N]
  cat(sprintf("\nrescued to conserved: %s confirmed intact at gametolog identity by miniprot",
              pc(k, nrow(r))))
  cat(sprintf("  |  positive control: %s  |  negative control: %.1f%%\n",
              pc(g[set == "conserved_Y" & mp == "intact, gametolog identity", .N], g[set == "conserved_Y", .N]),
              100 * fp))
  s <- d[cat_rescue == "degenerated"]
  cat(sprintf("still degenerated after the rescue: %s intact at gametolog identity (missed by TOGA2)\n",
              pc(s[mp == "intact, gametolog identity", .N], nrow(s))))
}
ii <- d[mp == "intact, gametolog identity"]
cat(sprintf("\n== the %d degenerated genes with an intact gametolog-identity hit ==\n", nrow(ii)))
cat(sprintf("  at the same neo-Y locus as a TOGA2 projection of the same Anax gene: %s\n",
            pc(ii[!is.na(toga_locus), .N], nrow(ii))))
cat(sprintf("  overlapping an annotated A. juncea neo-Y gene model:               %s\n",
            pc(ii[!is.na(y_model), .N], nrow(ii))))

# per-gene calls
setcolorder(g, c("set", "ag", "og", "jg", "mp"))
fwrite(g[, !c("ord", "ident0", "chr0")][order(set, mp, ag)], file.path(MP, "miniprot_calls.tsv"), sep = "\t")
cat("\nper-gene calls: ", file.path(MP, "miniprot_calls.tsv"), "\n")
