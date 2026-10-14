# One 11-tip CDS FASTA per locus of locus_table.tsv (per_locus/<pg>/<pg>.fasta); existing ones are skipped
suppressPackageStartupMessages({
  library(data.table)
  library(Biostrings)   # BiocManager::install("Biostrings") if missing
  library(ape)
})

locus_table_path <- "locus_table.tsv"
outdir           <- "per_locus"

# CDS FASTA per species; headers are gene IDs, with or without .t[0-9]+ (stripped).
# Aeshna_juncea provides both the neo-X and neo-Y tips
species_cds <- c(
  Aeshna_juncea          = "../divergence_dating/Aeshna.juncea.longest.cds.fna",
  Anax_junius            = "../divergence_dating/Anax.junius.longest.cds.fna",
  Somatochlora_arctica   = "../divergence_dating/Somatochlora.arctica.longest.cds.fna",
  Cordulegaster_boltonii = "../divergence_dating/Cordulegaster.boltonii.longest.cds.fna",
  Sympetrum_striolatum   = "../divergence_dating/Sympetrum.striolatum.longest.cds.fna",
  Orthetrum_coerulescens = "../divergence_dating/Orthetrum.coerulescens.longest.cds.fna",
  Pantala_flavescens     = "../divergence_dating/Pantala.flavescens.longest.cds.fna",
  Ischnura_elegans       = "../divergence_dating/Ischnura.elegans.longest.cds.fna",
  Platycnemis_pennipes   = "../divergence_dating/Platycnemis.pennipes.longest.cds.fna",
  Archilestes_grandis    = "../divergence_dating/Archilestes.grandis.longest.cds.fna"
)

gid <- function(x) sub("\\.t[0-9]+$", "", x)

locus_table <- fread(locus_table_path)
dir.create(outdir, showWarnings = FALSE)

# load each species' CDS, named by gene ID (stays NULL if the file is missing)
cat("=== Loading CDS FASTAs ===\n")
cds_by_species <- vector("list", length(species_cds))
names(cds_by_species) <- names(species_cds)
for (sp in names(species_cds)) {
  if (!file.exists(species_cds[[sp]])) next
  s <- readDNAStringSet(species_cds[[sp]])
  names(s) <- gid(sub("\\s.*", "", names(s)))
  cds_by_species[[sp]] <- s
}
missing_cds <- names(cds_by_species)[sapply(cds_by_species, is.null)]
if (length(missing_cds)) {
  cat("ERROR: missing/unset CDS paths for:", paste(missing_cds, collapse = ", "), "\n")
  cat("Edit `species_cds` at the top of the script.\n")
  quit(status = 1)
}
for (sp in names(cds_by_species)) {
  cat(sprintf("  %-30s %d sequences\n", sp, length(cds_by_species[[sp]])))
}

# loci still to extract
existing_fastas <- file.path(outdir, locus_table$pg_x,
                             paste0(locus_table$pg_x, ".fasta"))
to_extract <- !file.exists(existing_fastas)

if (any(to_extract)) {
  cat(sprintf("\n=== Stage 1: extracting %d per-locus FASTAs ===\n", sum(to_extract)))
  extraction_report <- data.table(pg_x = character(), status = character(),
                                  missing_tips = character())

  for (i in which(to_extract)) {
    row <- locus_table[i]
    pg  <- row$pg_x
    locus_dir <- file.path(outdir, pg)
    dir.create(locus_dir, showWarnings = FALSE)

    # tip label, species, gene ID
    tips <- list(
      list(tip = "Aeshna_juncea_neoX", sp = "Aeshna_juncea", g = row$g_x),
      list(tip = "Aeshna_juncea_neoY", sp = "Aeshna_juncea", g = row$g_y)
    )
    for (sp in setdiff(names(species_cds), "Aeshna_juncea")) {
      tips[[length(tips) + 1]] <- list(tip = sp, sp = sp, g = row[[sp]])
    }

    # sequence of each tip: exact gene ID first, then without the isoform suffix
    seqs <- list(); missing_tips <- character()
    for (t in tips) {
      cds <- cds_by_species[[t$sp]]
      hit <- cds[names(cds) == t$g]
      if (length(hit) == 0) hit <- cds[names(cds) == gid(t$g)]
      if (length(hit) == 0) missing_tips <- c(missing_tips, sprintf("%s(%s)", t$tip, t$g))
      else seqs[[t$tip]] <- as.character(hit[[1]])
    }

    # write the FASTA only if every tip was found
    if (length(missing_tips)) {
      extraction_report <- rbind(extraction_report,
        data.table(pg_x = pg, status = "incomplete",
                   missing_tips = paste(missing_tips, collapse = ",")))
    } else {
      ss <- DNAStringSet(unlist(seqs))
      names(ss) <- names(seqs)
      writeXStringSet(ss, file.path(locus_dir, paste0(pg, ".fasta")))
      extraction_report <- rbind(extraction_report,
        data.table(pg_x = pg, status = "complete", missing_tips = ""))
    }
  }
  fwrite(extraction_report, file.path(outdir, "extraction_report.tsv"), sep = "\t")

  n_complete <- sum(extraction_report$status == "complete")
  n_incomplete <- sum(extraction_report$status == "incomplete")
  cat(sprintf("Extracted: %d complete, %d incomplete\n", n_complete, n_incomplete))
  if (n_incomplete > 0) {
    cat("Incomplete loci (see extraction_report.tsv):\n")
    print(extraction_report[status == "incomplete"])
  }
} else {
  cat("\n=== Stage 1 skipped: all per-locus FASTAs already exist ===\n")
}
