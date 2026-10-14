# Per-gene depth, heterozygosity, callable bp and class -> genes.tsv (neo-X arm classes from gene_classes.tsv)
library(data.table)

XCHR <- "OZ369554.1"
BP   <- as.numeric(readLines("/path/to/project/gene_pairs_ds/breakpoint.txt"))   # fusion breakpoint; neo-X arm left, ancestral X right
NEOY <- "OZ369555.1"                  # 19.76x = 0.51x autosomal
EXCL <- "OZ369553.1"                  # 0.75x, unresolved - keep out of the autosomal baseline
MINLEN <- 10000000             # scaffolds shorter than this are treated as unplaced

fai <- fread("/path/to/project/data/genomes/Aeshna.juncea.clean.masked.fna.fai",
             select = 1:2, col.names = c("chr", "len"))
autosomes <- fai[len >= MINLEN & !chr %in% c(XCHR, NEOY, EXCL), chr]

# per-gene depth (length-weighted over CDS intervals), MAPQ 20 and MAPQ 0
cds <- fread("male.cds.regions.bed.gz", col.names = c("chr", "start", "end", "gene", "depth"))
cds[, len := end - start]
g <- cds[, .(chr = chr[1], start = min(start), cds_bp = sum(len),
             depth = sum(depth * len) / sum(len)), by = gene]
cds0 <- fread("male.cds.q0.regions.bed.gz",
              col.names = c("chr", "start", "end", "gene", "depth0"))
cds0[, len := end - start]
g0 <- cds0[, .(depth0 = sum(depth0 * len) / sum(len)), by = gene]
g <- merge(g, g0, by = "gene", all.x = TRUE)

# callable bp and het counts
g <- merge(g, fread("callable_cds_bp.tsv", col.names = c("gene", "callable_bp")), all.x = TRUE)
g <- merge(g, fread("het_per_gene.tsv",    col.names = c("gene", "het_n")),       all.x = TRUE)
g[is.na(callable_bp), callable_bp := 0]
g[is.na(het_n),       het_n       := 0]

# class: neo-X arm genes from gene_classes.tsv, autosomal and ancestral X from coordinates
g <- merge(g, fread("gene_classes.tsv", col.names = c("gene", "class")), all.x = TRUE)
g[chr %in% autosomes,          class := "autosomal"]
g[chr == XCHR & start >= BP,   class := "ancestral_X"]
g <- g[!is.na(class)]                       # drops neo-Y genes, unplaced, unclassified neo-X

# drop gene models rejected by TOGA2
dropped <- readLines("toga2_dropped_genes.txt")
cat(sprintf("[toga2] %d dropped models found in the table, removed\n", sum(g$gene %chin% dropped)))
g <- g[!gene %chin% dropped]

# depth relative to the autosomal median
g[, norm_depth  := depth  / median(depth[class == "autosomal"])]
g[, norm_depth0 := depth0 / median(depth0[class == "autosomal"])]
g[, mapq_loss   := 1 - depth / depth0]

# het sites per callable kb
g[, het_kb := het_n / (callable_bp / 1000)]
g[callable_bp < 300, het_kb := NA]          # too little callable CDS to score

print(g[, .(n = .N, med_depth = round(median(norm_depth), 3),
            med_het = round(median(het_kb, na.rm = TRUE), 2)), by = class])
fwrite(g, "genes.tsv", sep = "\t")