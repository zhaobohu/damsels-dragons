# Name fixes the GBIF resolver (00_resolve_names.R) cannot make, and a suggester for unplaced names
suppressPackageStartupMessages(library(data.table))

# OCR-damaged karyotype names GBIF cannot match, repaired before the query
PRE_FIX <- c(
  "Brachvmesia furcata" = "Brachymesia furcata",
  "Brachvmesia gravida" = "Brachymesia gravida",
  "Brachvmesia herbida" = "Brachymesia herbida"
)

# GBIF accepted name -> name to use, applied after resolution; isoceles is moved
# out of Aeshna to join the backbone tip Isoaeshna isoceles
GBIF_OVERRIDE <- c(
  "Aeshna isoceles" = "Isoaeshna isoceles"
)

# family misspellings in the karyotype table
FAMILY_FIX <- c(
  "Gordulegastridae" = "Cordulegastridae",
  "Platysticitidae"  = "Platystictidae"
)

# genus, epithet, and indeterminate-name test
genus_of <- function(x) sub(" .*", "", x)
epi_of   <- function(x) sub("^\\S+ ", "", x)
is_indet <- function(x) grepl("\\bsp\\.?$|\\bspp\\.?$|\\bnr\\b", x, ignore.case = TRUE)

# suggest resolved names (pool) for unresolved ones and write them to name_review.tsv
# for review by hand; nothing is applied automatically
suggest_names <- function(unresolved, un_fam, pool, pool_fam,
                          review = "name_review.tsv") {
  un <- unique(data.table(query = unresolved, fam = un_fam)[!is_indet(query)])
  pl <- unique(data.table(sp = pool, fam = pool_fam)[!is_indet(sp)])
  sug <- list()

  # within-genus spelling variants (epithet edit distance 1-2)
  for (g in unique(genus_of(un$query))) {
    mm <- un$query[genus_of(un$query) == g]
    pp <- pl$sp[genus_of(pl$sp) == g]
    if (!length(pp)) next
    d <- adist(epi_of(mm), epi_of(pp))
    for (r in seq_along(mm)) {
      j <- which(d[r, ] <= 2 & d[r, ] > 0)
      if (length(j)) sug[[length(sug) + 1]] <-
        data.table(query = mm[r], suggested = pp[j], method = "spelling", dist = d[r, j])
    }
  }

  # same epithet and family, different genus (may pair different species)
  j2 <- merge(data.table(query = un$query, fam = un$fam, ep = epi_of(un$query)),
              data.table(suggested = pl$sp, fam = pl$fam, ep = epi_of(pl$sp)),
              by = c("fam", "ep"), allow.cartesian = TRUE)
  j2 <- j2[genus_of(query) != genus_of(suggested)]
  if (nrow(j2)) sug[[length(sug) + 1]] <-
    j2[, .(query, suggested, method = "reassignment", dist = NA_integer_)]

  # keep earlier accept decisions and write the review table
  if (length(sug)) {
    sg <- unique(rbindlist(sug))[order(query, method, dist)]
    sg[, accept := 0L]
    if (file.exists(review)) {
      old <- fread(review)
      if ("accept" %in% names(old)) sg[old, on = .(query, suggested), accept := i.accept]
    }
    fwrite(sg, review, sep = "\t")
    message(sprintf("  %d suggestions in %s", nrow(sg), review))
  }
  invisible(NULL)
}
