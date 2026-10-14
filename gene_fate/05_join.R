# TOGA2 statuses per Anax gene (by where each copy lands: neo-Y / X / elsewhere) and per A. juncea gene
# (valid / broken lists for the pangenes + MCScanX-pairs filter); writes to <run dir>/toga2_out/join/
# usage: Rscript 05_join.R <run dir> <A. juncea gene GFF>
suppressMessages(library(data.table))
a <- commandArgs(TRUE); RUN <- a[1]; QRY_GFF <- a[2]
OUT   <- file.path(RUN, "toga2_out"); REFIN <- file.path(RUN, "ref_input")
CHAIN <- Sys.glob(file.path(RUN, "chains/results/07_final/*.chain.gz"))[1]
REN   <- file.path(RUN, "genomes/Aeshna_juncea.rename.tsv")
RAWBED <- Sys.glob(file.path(OUT, "meta/query_annotation.with_discarded_exons.bed*"))[1]  # all aligned projections
X <- "OZ369554"; Y <- "OZ369555"          # A. juncea fused X and neo-Y, as named in the chains (no ".1")
BRK <- as.numeric(readLines("/path/to/project/gene_pairs_ds/breakpoint.txt"))   # positions < BRK = neo-X arm
dir.create(file.path(OUT, "join"), showWarnings = FALSE)

# TOGA2's ranking (cesar_wrapper_constants.py): FI > I > PI > UL > L > M > PM > PG > PP > N
rk   <- c(N = -1, PP = 0, PG = 1, PM = 2, M = 3, L = 4, UL = 5, PI = 6, I = 7, FI = 8)
best <- function(s) { i <- which.max(rk[s]); if (length(i)) s[i] else NA_character_ }

# projections (all aligned ones, not only the final annotation): <transcript>#<chain>[,<chain>...],
# transcript = <ref transcript>#<ref gene>
lss <- fread(file.path(OUT, "meta/loss_summary_extended.tsv"))         # level  entry  status
pr  <- lss[level == "PROJECTION", .(proj = entry, status)]
pr[, base   := sub("#(paralog|retro)$", "", sub("\\$.*$", "", proj))]  # drop fragment / label suffixes
pr[, chains := sub(".*#", "", base)]
pr[, tr     := sub("#[^#]*$", "", base)]

# paralog-labelled copies come back as PG instead of L / M: undo that with TOGA2's own rules
# (transcript_meta.tsv: intact < 0.35, non-deleted < 0.2 or middle affected -> L, else M)
TM <- Sys.glob(file.path(OUT, "meta/transcript_meta.tsv*"))[1]
tm <- fread(cmd = sprintf("zcat -f %s", shQuote(TM)), select = c(1, 6, 7, 8),
            col.names = c("proj", "nondel", "intact", "middle"))
pr <- merge(pr, tm, by = "proj", all.x = TRUE)
pr[, was_pg := status == "PG"]
pr[was_pg & !is.na(middle),
   status := fifelse(intact < 0.35 | nondel < 0.2 | middle == "MIDDLE_IS_AFFECTED", "L", "M")]
cat(sprintf("[PG] %d paralog-labelled projections re-assigned: %d L, %d M; %d left as PG (no metadata)\n",
            pr[was_pg & status != "PG", .N], pr[was_pg & status == "L", .N],
            pr[was_pg & status == "M", .N], pr[status == "PG", .N]))
iso <- fread(file.path(REFIN, "toga.isoforms.tsv"), header = FALSE, col.names = c("gene", "tr"))
rb  <- fread(file.path(REFIN, "toga.transcripts.bed"), header = FALSE, select = c(1, 4),
             col.names = c("ref_chr", "tr"))

# where each projection lands: chain id -> query sequence
ch <- fread(cmd = sprintf("zcat %s | awk '$1==\"chain\"{print $13\"\\t\"$8}'", shQuote(CHAIN)),
            header = FALSE, col.names = c("chain", "qchr"), colClasses = "character")
sp <- strsplit(pr$chains, ",", fixed = TRUE)                           # fragmented: several chains
pc <- pr[rep(seq_len(.N), lengths(sp)), .(proj, base, tr, status)][, chain := unlist(sp)]
pc <- merge(pc, ch, by = "chain", all.x = TRUE)
if (pc[is.na(qchr), .N]) warning(pc[is.na(qchr), .N], " projection chain ids not found in ", CHAIN)
reg <- pc[, .(onY = any(qchr == Y, na.rm = TRUE), onX = any(qchr == X, na.rm = TRUE)),
          by = .(proj, base, tr, status)]
reg[, region := fifelse(onY, "Y", fifelse(onX, "X", "other"))]
reg <- merge(reg, iso, by = "tr")

# (a) per Anax gene: best status, number and all statuses of its copies on the neo-Y, X, elsewhere
per <- reg[, .(best = best(status), n = uniqueN(base), all = paste(sort(unique(status)), collapse = ",")),
           by = .(gene, region)]
g <- merge(unique(iso[, .(gene)]),
           merge(iso, rb, by = "tr")[, .(ref_chr = ref_chr[1]), by = gene], by = "gene", all.x = TRUE)
g <- merge(g, lss[level == "GENE", .(gene = entry, genome_wide = status)], by = "gene", all.x = TRUE)
for (r in c("Y", "X", "other")) {
  x <- per[region == r, .(gene, best, n, all)]
  setnames(x, c("best", "n", "all"), paste0(r, c("_best", "_n", "_all")))
  g <- merge(g, x, by = "gene", all.x = TRUE)
}
g[is.na(Y_n), Y_n := 0L][is.na(X_n), X_n := 0L][is.na(other_n), other_n := 0L]
callf <- function(n, b, none) fcase(n == 0,                  none,
                                    b %chin% c("FI", "I"),   "intact",
                                    b == "PI",               "partially_intact",
                                    b %chin% c("UL", "L"),   "inactivated",
                                    b %chin% c("M", "PM"),   "missing",
                                    default = "other")
g[, `:=`(Y_call = callf(Y_n, Y_best, "no_Y_copy"), X_call = callf(X_n, X_best, "no_X_copy"))]
fwrite(g, file.path(OUT, "join/anax_gene_status.tsv"), sep = "\t")

# Anax genes with no intact copy anywhere: a lost copy can't be told from a bad model
ok  <- per[best %chin% c("FI", "I", "PI"), unique(gene)]
unv <- g[!gene %chin% ok]
writeLines(unv$gene, file.path(OUT, "join/anax_unverified_genes.txt"))
cat("\n== Anax genes with no intact copy anywhere (unverified), by Anax chromosome group ==\n")
print(g[, .(genes = .N, unverified = sum(!gene %chin% ok), pct = round(100 * mean(!gene %chin% ok), 1)),
        by = .(grp = fcase(ref_chr == "OZ349811", "neoX_homolog", ref_chr == "OZ349813", "ancX_homolog",
                           grepl("^OZ", ref_chr), "autosome", default = "scaffold"))])
# the X copy of the same gene is the matched baseline: how often TOGA2 calls a retained copy lost
cat("\n== genes on the Anax neo-X homolog (OZ349811): X copy (rows) x neo-Y copy (columns) ==\n")
print(g[ref_chr == "OZ349811", table(X_call, Y_call)])

# (b) per A. juncea gene. Per projection: best status over fragments, paralog label (L/M came back
# as PG) and the Anax gene it projects
pst <- merge(pr, iso, by = "tr")[, .(status = best(status), pg = any(was_pg), ag = gene[1]), by = base]
pst <- merge(pst, reg[, .(region = region[1]), by = base], by = "base", all.x = TRUE)
# gametolog rule: for a neo-X-homolog gene intact on one sex chromosome, a paralog-labelled L copy
# on the other one counts as evidence; paralog-labelled L is ignored everywhere else
gamY <- g[ref_chr == "OZ349811" & X_call == "intact", gene]      # X intact -> Y copy is the gametolog
gamX <- g[ref_chr == "OZ349811" & Y_call == "intact", gene]
pst[, gam := pg & ((region == "Y" & ag %chin% gamY) | (region == "X" & ag %chin% gamX))]

# projected exons (BED12 blocks) on juncea coordinates; overlap = juncea CDS on the same strand
# (gene spans are not enough: long introns span unrelated genes)
bed <- fread(cmd = sprintf("zcat -f %s | cut -f1-4,6,11,12", shQuote(RAWBED)), header = FALSE,
             col.names = c("chr", "s0", "e0", "proj", "strand", "sizes", "starts"))
bed[, base := sub("#(paralog|retro)$", "", sub("\\$.*$", "", proj))]
bed <- merge(bed, pst, by = "base")
if (!nrow(bed)) stop("no projection in ", RAWBED, " matched meta/loss_summary_extended.tsv")
bs <- strsplit(bed$starts, ",", fixed = TRUE); bz <- strsplit(bed$sizes, ",", fixed = TRUE)
bed[, nb := lengths(bs)]                                   # projected exons in this BED row
ex <- bed[rep(seq_len(.N), lengths(bs)), .(chr, strand, base, status, pg, ag, gam, s0)]
ex[, `:=`(start = s0 + as.integer(unlist(bs)) + 1L,
          end = s0 + as.integer(unlist(bs)) + as.integer(unlist(bz)))]

# juncea genes and CDS from the GFF, renamed as in the chains
gff <- fread(QRY_GFF, header = FALSE, sep = "\t", quote = "")
tag <- function(x, k) sub(sprintf("^(.*;)?%s=([^;]+).*$", k), "\\2", x)
ren <- fread(REN, header = FALSE, col.names = c("old", "new"))
gff[, V1 := ren$new[match(V1, ren$old)]]                      # OZ369555.1 (GFF) -> OZ369555 (TOGA2)
gf  <- gff[V3 == "gene", .(chr = V1, start = V4, end = V5, gene = tag(V9, "ID"))]
tx  <- gff[V3 %chin% c("mRNA", "transcript"), .(tx = tag(V9, "ID"), gene = tag(V9, "Parent"))]
cds <- merge(gff[V3 == "CDS", .(chr = V1, start = V4, end = V5, strand = V7, tx = tag(V9, "Parent"))],
             tx, by = "tx")
if (!nrow(cds)) stop("no CDS line in ", QRY_GFF, " links to a gene through its mRNA's Parent=")
setkey(cds, chr, start, end); setkey(ex, chr, start, end)
ov <- foverlaps(ex, cds, nomatch = NULL)[strand == i.strand, .(gene, chr, base, status, pg, ag, gam)]

# anax_projections.tsv: one row per (projection, query sequence) with Anax gene, region, status,
# paralog / retro labels, projected vs reference exons, overlapping juncea genes, coordinates
rex <- fread(file.path(REFIN, "toga.transcripts.bed"), header = FALSE, select = c(4, 10),
             col.names = c("tr", "ref_exons"))
pa <- bed[, .(start = min(s0) + 1L, end = max(e0), exons = sum(nb)), by = .(base, chr)]
pa <- merge(pa, bed[, .(para = any(grepl("#paralog$", proj)), retro = any(grepl("#retro$", proj))),
                    by = base], by = "base")
pa <- merge(pa, pst[, .(base, ag, status, pg)], by = "base")
pa <- merge(pa, merge(unique(pr[, .(base, tr)]), rex, by = "tr")[, .(ref_exons = max(ref_exons)), by = base],
            by = "base", all.x = TRUE)
pa <- merge(pa, ov[, .(juncea = paste(sort(unique(gene)), collapse = ",")), by = .(base, chr)],
            by = c("base", "chr"), all.x = TRUE)
pa[, para := para | pg][, pg := NULL]                      # L/M paralogs lost the label to PG; restore it
pa[, region := fcase(chr == Y, "Y", chr == X & start < BRK, "neoX", chr == X, "ancX",
                     grepl("^OZ", chr), "auto", default = "other")]
setcolorder(pa, c("ag", "region", "status", "para", "retro", "exons", "ref_exons", "juncea",
                  "chr", "start", "end", "base"))
fwrite(pa[order(ag, region)], file.path(OUT, "join/anax_projections.tsv"), sep = "\t")

# per juncea gene: broken only if nothing overlapping it beats L and an L comes from an orthologous
# projection of a verified Anax gene (l_src) or from a gametolog (gam_src);
# kept_reason for other L-best genes: paralog_L, unverified_L or both
js <- ov[, .(best = best(status), n_proj = uniqueN(base),
             l_src   = paste(sort(unique(ag[status == "L" & !pg & ag %chin% ok])), collapse = ","),
             gam_src = paste(sort(unique(ag[status == "L" & gam])), collapse = ","),
             pg_l  = any(status == "L" & pg & !gam), unv_l = any(status == "L" & !pg & !ag %chin% ok)),
         by = gene]
js[l_src == "", l_src := NA][gam_src == "", gam_src := NA]
js <- merge(gf, js, by = "gene", all.x = TRUE)
js[, grp := fcase(chr == Y, "neo-Y", chr == X & start < BRK, "X_neoX-arm", chr == X, "X_ancX-arm",
                  grepl("^OZ", chr), chr, default = "scaffolds")]
# the X side of the gametolog rule is for the neo-X arm only: on the anc-X arm it is an ordinary paralog
js[grp == "X_ancX-arm" & !is.na(gam_src), `:=`(gam_src = NA, pg_l = TRUE)]
js[, l_evidence := !is.na(l_src) | !is.na(gam_src)]
js[best == "L" & !l_evidence, kept_reason := fcase(pg_l & !unv_l, "paralog_L", unv_l & !pg_l, "unverified_L",
                                                   default = "both")]
js[, c("pg_l", "unv_l") := NULL]
js[is.na(best), `:=`(best = "no_projection", n_proj = 0L, l_evidence = FALSE)]
cat(sprintf("\n[gametolog] %d juncea genes dropped on gametolog evidence alone (%d on the neo-Y, %d on the neo-X arm)\n",
            js[best == "L" & is.na(l_src) & !is.na(gam_src), .N],
            js[best == "L" & is.na(l_src) & !is.na(gam_src) & chr == Y, .N],
            js[best == "L" & is.na(l_src) & !is.na(gam_src) & grp == "X_neoX-arm", .N]))
setcolorder(js, c("gene", "chr", "start", "end", "grp"))
setorder(js, chr, start)
fwrite(js, file.path(OUT, "join/juncea_gene_status.tsv"), sep = "\t")
writeLines(js[best %chin% c("FI", "I", "PI"), gene], file.path(OUT, "join/juncea_valid_genes.txt"))
writeLines(js[best == "L" & l_evidence,      gene], file.path(OUT, "join/juncea_broken_genes.txt"))

# summaries
cat("\n== juncea genes by best overlapping projection ==\n")
print(js[, .N, by = best][order(-N)])
cat(sprintf("\n[broken] %d of %d L-best juncea genes dropped; %d kept (L only from paralog-labelled\n",
            js[best == "L" & l_evidence, .N], js[best == "L", .N], js[best == "L" & !l_evidence, .N]),
    "         projections or from Anax genes with no intact copy anywhere)\n", sep = "")
cat("\n== L-best genes kept, by reason and region ==\n")
print(dcast(js[best == "L" & !l_evidence, .N,
               by = .(grp = fifelse(grp %chin% c("neo-Y", "X_neoX-arm", "X_ancX-arm"), grp, "elsewhere"),
                      kept_reason)],
            grp ~ kept_reason, value.var = "N", fill = 0L))
cat("\n== juncea genes by chromosome (unplaced scaffolds pooled) ==\n")
print(js[, .(genes = .N, pct_valid = round(100 * mean(best %chin% c("FI", "I", "PI")), 1),
             L_best = sum(best == "L"), dropped = sum(best == "L" & l_evidence)), by = grp][
  , pct_dropped := round(100 * dropped / genes, 1)][order(grp)])
# runs of >= 3 adjacent dropped genes: few distinct Anax genes behind a run = one family / tandem array
cl <- js[, .(chr, start, end, l_src, gam_src, drop = best == "L" & l_evidence)][, run := rleid(chr, drop)][
  drop == TRUE, .(chr = chr[1], start = min(start), end = max(end), dropped = .N,
                  anax_genes = uniqueN(na.omit(unlist(strsplit(c(l_src, gam_src), ","))))),
  by = run][dropped >= 3][order(-dropped)]
cat(sprintf("\n== %d runs of >= 3 adjacent dropped genes (%d genes) ==\n", nrow(cl), sum(cl$dropped)))
print(cl[, !"run"])
