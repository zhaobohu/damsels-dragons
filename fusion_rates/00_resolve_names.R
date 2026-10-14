#!/usr/bin/env Rscript
# Resolve karyotype, BOLD and backbone names against the GBIF backbone -> name_map.csv

suppressPackageStartupMessages({ library(data.table); library(rgbif); library(ape) })
source("config.R")
source("reconcile.R")     # PRE_FIX, FAMILY_FIX, GBIF_OVERRIDE

# karyotype table (family spellings fixed) and BOLD records
k <- fread(KARYO_CSV)
ff <- FAMILY_FIX
k[tax3 %chin% names(ff), tax3 := ff[tax3]]

b <- fread(BOLD_TSV, quote = "", sep = "\t")

# backbone tips are ID_Family_Genus_species
bb  <- read.nexus(BACKBONE_TRE)
pb  <- tstrsplit(gsub("^'|'$", "", bb$tip.label), "_", fixed = TRUE)
bbt <- data.table(bb_tip = bb$tip.label, family = pb[[2]],
                  genus = pb[[3]], species = paste(pb[[3]], pb[[4]]))
fwrite(bbt, "backbone_families.tsv", sep = "\t")

# every name, tagged by source
raw <- unique(rbindlist(list(
  data.table(source = "karyo",    verbatim = k$sp),
  data.table(source = "bold",     verbatim = b$species),
  data.table(source = "backbone", verbatim = bbt$species)
))[nzchar(trimws(verbatim))])

# repair OCR-damaged strings before the query
raw[, query := ifelse(verbatim %chin% names(PRE_FIX), PRE_FIX[verbatim], verbatim)]

u <- unique(raw$query)
cat(sprintf("%d unique names (%d karyotype, %d BOLD, %d backbone)\n", length(u),
            uniqueN(raw[source == "karyo", query]),
            uniqueN(raw[source == "bold", query]),
            uniqueN(raw[source == "backbone", query])))

# GBIF match; stop with a clear error if expected rgbif columns are missing
g <- name_backbone_checklist(u, kingdom = "Animalia")
stopifnot(nrow(g) == length(u))

need <- c("verbatim_name", "species", "speciesKey", "genus", "family",
          "order", "status", "matchType", "confidence")
miss <- setdiff(need, names(g))
if (length(miss))
  stop("rgbif columns missing: ", paste(miss, collapse = ", "),
       "\navailable: ", paste(names(g), collapse = ", "))

# keep species-level matches only
gm <- as.data.table(g)[, .(query = verbatim_name, acc = species, gbif_key = speciesKey,
                           gen = genus, fam = family, ord = order, st = status,
                           mt = matchType, conf = confidence)]
gm[mt %chin% c("HIGHERRANK", "NONE") | is.na(gbif_key), `:=`(acc = NA_character_, gbif_key = NA)]

# reject matches outside ALLOWED_ORDERS (e.g. Cordulia amurensis matched to a bivalve)
ok_ord <- ALLOWED_ORDERS
bad <- gm[!is.na(gbif_key) & (is.na(ord) | !ord %chin% ok_ord)]
if (nrow(bad)) {
  cat("\nREJECTED -- GBIF matched these outside ",
      paste(ok_ord, collapse = "/"), ":\n", sep = "")
  print(bad[, .(query, acc, gen, fam, ord, st, conf)])
  gm[!is.na(gbif_key) & (is.na(ord) | !ord %chin% ok_ord),
     `:=`(acc = NA_character_, gbif_key = NA)]
}

# apply GBIF_OVERRIDE, then give each final name one gbif_key so overridden names join their target
ov <- GBIF_OVERRIDE
n_ov <- gm[!is.na(acc) & acc %chin% names(ov), .N]
gm[!is.na(acc) & acc %chin% names(ov), acc := ov[acc]]
gm[!is.na(acc), gbif_key := min(gbif_key, na.rm = TRUE), by = acc]
if (n_ov) cat(sprintf("\n%d names remapped by GBIF_OVERRIDE\n", n_ov))

m <- merge(raw, gm, by = "query", all.x = TRUE)
setcolorder(m, c("source", "verbatim", "query", "acc", "gbif_key"))
fwrite(m, "name_map.csv")

# report: match types, unresolved, fuzzy matches, synonyms and genus reassignments
cat("\nmatchType x status:\n"); print(table(m$mt, m$st, useNA = "ifany"))

cat(sprintf("\nunresolved: %d of %d\n", m[is.na(acc), .N], nrow(m)))
print(m[is.na(acc), .N, by = source])

cat("\nFUZZY matches -- INSPECT THESE BY HAND:\n")
print(unique(m[mt == "FUZZY", .(source, verbatim, acc, conf, st)])[order(conf)])

cat("\nnames collapsing onto a shared accepted name (synonyms resolved):\n")
cl <- unique(m[!is.na(acc), .(verbatim, acc, st, source)])
print(cl[acc %chin% cl[, .N, by = acc][N > 1, acc]][order(acc)])

cat("\ngenus reassignments (verbatim genus != GBIF genus):\n")
print(unique(m[!is.na(acc) & sub(" .*", "", verbatim) != gen,
               .(verbatim, acc, st)])[order(verbatim)])

# karyotype species that share a GBIF key with BOLD or the backbone
kk <- m[source == "karyo" & !is.na(gbif_key), unique(gbif_key)]
av <- m[source %chin% c("bold", "backbone") & !is.na(gbif_key), unique(gbif_key)]
cat(sprintf("\nkaryotype species with a GBIF gbif_key: %d\n", length(kk)))
cat(sprintf("  ...of which also present in BOLD or the backbone: %d\n",
            sum(kk %in% av)))

# names moving into or out of the focal genus (these change the clade under test)
fg <- FOCAL_GENUS
moved <- unique(m[!is.na(acc) &
                  (sub(" .*", "", verbatim) == fg) != (sub(" .*", "", acc) == fg),
                  .(source, verbatim, acc, st)])
cat(sprintf("\n=== names moving INTO or OUT OF %s -- decide these deliberately ===\n", fg))
if (nrow(moved)) {
  print(moved[order(verbatim)])
  cat("Add to GBIF_OVERRIDE in reconcile.R to keep the source name.\n")
} else cat("none\n")

# karyotype names GBIF could not resolve
cat("\nkaryotype names GBIF could not resolve (add to PRE_FIX in reconcile.R\n")
cat("if they are typos, otherwise they are genuinely unplaceable):\n")
print(m[source == "karyo" & is.na(acc), verbatim])
