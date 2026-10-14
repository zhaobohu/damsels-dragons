#!/usr/bin/env Rscript
# Join karyotypes, BOLD COI and the backbone on GBIF keys (name_map.csv) -> pipeline inputs

suppressPackageStartupMessages({ library(data.table); library(ape) })
source("config.R")
source("reconcile.R")

if (!file.exists("name_map.csv"))
  stop("name_map.csv not found -- run 00_resolve_names.R first")
nm <- fread("name_map.csv")

# verbatim name -> gbif_key per source, and gbif_key -> accepted name
# (the column is gbif_key: data.table(key = ...) would be taken as the setkey argument)
key_of <- function(src, v) nm[source == src][match(v, verbatim), gbif_key]
acc_of <- function(k)      nm[!is.na(gbif_key)][match(k, gbif_key), acc]

# karyotypes, family spellings fixed
k <- fread(KARYO_CSV)
k <- k[nzchar(trimws(sp))]
ff <- FAMILY_FIX
k[tax3 %chin% names(ff), tax3 := ff[tax3]]

# drop records with an empty sex field (counted for the report)
n_blank  <- k[!nzchar(trimws(sex)), .N]
sp_blank <- k[, .(ab = all(!nzchar(trimws(sex)))), by = sp][ab == TRUE, .N]
k <- k[nzchar(trimws(sex))]

# rename to the GBIF accepted name; drop names GBIF could not place
k[, sp_orig := sp]
k[, gbif_key := key_of("karyo", sp_orig)]
k[!is.na(gbif_key), sp := acc_of(gbif_key)]
n_renamed <- k[sp != sp_orig, uniqueN(sp_orig)]
k <- k[!is.na(gbif_key)]
fwrite(unique(k[sp != sp_orig, .(sp_orig, sp)])[order(sp_orig)],
       "name_applied.tsv", sep = "\t")

# haploid linkage groups = diploid autosomes (aut) / 2 + X elements (two for X1X2Y)
k[, `:=`(neo = sex %chin% c("neo-XY", "X1X2Y", "XY"),
         nchr_rec = aut / 2 + fifelse(sex == "X1X2Y", 2L, 1L))]

# one row per species: any neo record makes it neo-sex, and nchr is the mode among those records
sp <- k[, {
  use <- if (any(neo)) .SD[neo == TRUE] else .SD
  m   <- as.numeric(names(sort(table(use$nchr_rec), decreasing = TRUE)))[1]
  .(gbif_key = gbif_key[1], family = tax3[1], suborder = con[1],
    sex_system = if (any(neo)) sort(unique(use$sex))[1] else "X",
    nchr = m, n_rec = .N,
    conflict = uniqueN(na.omit(nchr_rec)) > 1 || uniqueN(neo) > 1)
}, by = .(species = sp)]
sp[, genus := sub(" .*", "", species)]

# records of species with disagreeing karyotypes
cf <- sp[conflict == TRUE, species]
fwrite(k[sp %chin% cf, .(sp, aut, sex, loc, karyo.source)][order(sp)],
       "karyo_conflicts.tsv", sep = "\t")

# backbone
bb  <- read.nexus(BACKBONE_TRE)
bbt <- fread("backbone_families.tsv")     # written by 00_resolve_names.R
bbt[, gbif_key := key_of("backbone", species)]

# BOLD COI-5P records of karyotyped species, kept at record level
b <- fread(BOLD_TSV, quote = "", sep = "\t")
b <- b[marker_code == "COI-5P" & nchar(nuc) > 0 & nzchar(species)]
b[, gbif_key := key_of("bold", species)]
b <- b[gbif_key %in% sp$gbif_key]
b[, `:=`(nuc = gsub("^-*(.*?)-*$", "\\1", nuc),
         record_id = gsub("[^A-Za-z0-9]+", "_", processid),
         acc_name = acc_of(gbif_key))]
stopifnot(!anyDuplicated(b$record_id))

# write the records and their metadata
s <- character(0)
for (i in seq_len(nrow(b)))
  s <- c(s, sprintf(">%s", b$record_id[i]), gsub("-", "", b$nuc[i]))
writeLines(s, "unfiltered.COI.fna")
fwrite(b[, .(record_id, gbif_key, acc_name, species_bold = species,
             family_bold = family, bin_uri, insdc_acs, source = "bold")],
       "prep.csv")

# BIN discordance: species in several BINs, BINs shared by several species
bd <- b[nzchar(bin_uri), .(n_bin = uniqueN(bin_uri)), by = acc_name][n_bin > 1]
bs <- b[nzchar(bin_uri), .(n_sp = uniqueN(acc_name),
                           spp = paste(sort(unique(acc_name)), collapse = "; ")),
        by = bin_uri][n_sp > 1]
fwrite(bd, "bin_flags.tsv", sep = "\t"); fwrite(bs, "bin_shared.tsv", sep = "\t")

# species still needing COI from NCBI, and the taxon table
sp[, has_coi := gbif_key %in% b$gbif_key]
writeLines(sp[has_coi == FALSE, species], "needed_from_ncbi.txt")
fwrite(sp[, .(species, genus, family, nchr, sex_system, suborder, gbif_key, has_coi, conflict)],
       TAXA_TSV, sep = "\t")

# suggestions for karyotype names GBIF could not place (pool and families aligned row-for-row)
kraw <- fread(KARYO_CSV)[nzchar(trimws(sp))]
kraw[tax3 %chin% names(ff), tax3 := ff[tax3]]
unres <- unique(kraw[is.na(key_of("karyo", sp)), .(sp, tax3)])
if (nrow(unres)) {
  pool <- unique(rbind(b[, .(sp = acc_name, fam = family)],
                       bbt[!is.na(gbif_key), .(sp = species, fam = family)]))
  invisible(suggest_names(unres$sp, unres$tax3, pool$sp, pool$fam))
}

# report
cat("=== karyotypes ===\n")
cat(sprintf("dropped %d records with an empty sex field (%d species lost entirely)\n",
            n_blank, sp_blank))
cat(sprintf("%d species renamed to their GBIF accepted name\n", n_renamed))
cat(sprintf("%d records -> %d species (%d with conflicting records)\n",
            nrow(k), nrow(sp), sp[, sum(conflict)]))
print(sp[, .(species = .N, neo = sum(sex_system != "X")), by = suborder])

cat("\n=== COI coverage by family ===\n")
cv <- sp[, .(karyo = .N, coi = sum(has_coi)), by = .(suborder, family)][order(-karyo)]
cv[, pct := round(100 * coi / karyo)]
print(cv[karyo >= 3], nrows = 30)
cat(sprintf("\noverall: %d of %d karyotyped species have BOLD COI (%.0f%%)\n",
            sp[, sum(has_coi)], nrow(sp), 100 * sp[, mean(has_coi)]))
cat(sprintf("%d BOLD records kept for those species\n", nrow(b)))

cat("\n=== the number that limits the analysis ===\n")
cat(sprintf("neo-sex species in karyotype data : %d\n", sp[sex_system != "X", .N]))
cat(sprintf("  ...of those, with COI            : %d\n", sp[sex_system != "X" & has_coi, .N]))
cat(sprintf("  ...of those, in Aeshna           : %d\n",
            sp[sex_system != "X" & has_coi & genus == FOCAL_GENUS, .N]))

cat("\n=== backbone ===\n")
cat(sprintf("%d tips, %d families, root depth %.2f time units\n",
            Ntip(bb), uniqueN(bbt$family), max(branching.times(bb))))
cat(sprintf("Aeshnidae tips: %d.  %s tips: %d\n", bbt[family == "Aeshnidae", .N],
            FOCAL_GENUS, bbt[genus == FOCAL_GENUS, .N]))
cat(sprintf("karyotyped species also on the backbone (by GBIF key): %d of %d\n",
            sp[gbif_key %in% bbt$gbif_key, .N], nrow(sp)))
