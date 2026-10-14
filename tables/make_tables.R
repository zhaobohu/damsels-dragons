#!/usr/bin/env Rscript
# Supplementary Tables S1-S5 (TSV, XLSX, table_notes.tsv); runs in tables/, inputs relative to the repository root

suppressPackageStartupMessages({ library(data.table); library(ape) })

RUNS <- ".."                  # the repository root
FIG4 <- "../figs/fig4/fig4.R" # re-sourced for Table S4

# the nine de-novo annotated genomes: BRAKER3/GALBA script code -> genome (= run directory under annotation/genes);
# the tenth, Ischnura elegans, has the NCBI RefSeq annotation
CODES <- c(AeshJunc = "Aeshna_juncea", AnaxJuni = "Anax_junius", ArchGran = "Archilestes_grandis",
           CordBolt = "Cordulegaster_boltonii", OrthCoer = "Orthetrum_coerulescens",
           PantFlav = "Pantala_flavescens", PlatPenn = "Platycnemis_pennipes",
           SomaArct = "Somatochlora_arctica", SympStri = "Sympetrum_striolatum")
BUSCO_TXT <- "short_summary.specific.insecta_odb12.2.busco_combined.txt"

# inputs, relative to RUNS
IN <- list(
  assemblies   = "palaeoptera_synteny/pairwise_synteny/assemblies.tsv", # the 29 assemblies of Fig. 1
  kmer         = "microchromosome/sex_kmers/kmer.sh",                  # Ae. juncea read and female-assembly accessions
  annot_sbatch = sprintf("annotation/genes/%s/%s.gene.annot.sbatch", CODES, names(CODES)),
  braker_busco = sprintf("annotation/genes/%s/busco_combined/%s", CODES, BUSCO_TXT),
  annevo_busco = "palaeoptera_synteny/genespace_run/annevo.busco.csv",
  annevo_genes = "palaeoptera_synteny/annevo_annotation/gene.counts.tsv",
  mcmc_ctl     = "divergence_dating/full_chain1/mcmctree.ctl",         # the settings of the full-data run
  calib_tree   = "divergence_dating/ajun_calib.tre",
  ages         = "divergence_dating/modelsel/ages_all_runs.tsv",       # incl. the prior-only run
  n0           = "palaeoptera_synteny/genespace_run/genespace/orthofinder/Results_Aug18/Phylogenetic_Hierarchical_Orthogroups/N0.tsv",
  micro_beds   = sprintf("palaeoptera_synteny/genespace_run/genespace/bed/%s.bed",
                         c("Ischnura_elegans", "Platycnemis_pennipes", "Tanypteryx_hageni")),
  juncea_cds   = "microchromosome/cds_depth/cds.annevo.bed")
p <- function(f) file.path(RUNS, f)

# TRUE if every named input exists, else names the missing ones (the table is skipped)
have <- function(...) {
  f <- unlist(IN[c(...)]); m <- f[!file.exists(p(f))]
  if (length(m)) message("  missing (skipped): ", paste(file.path(RUNS, m), collapse = ", "))
  length(m) == 0
}
# compare a table value with the manuscript text: DIFF means one of them needs changing
CHECK <- function(label, got, want) {
  ok <- isTRUE(all.equal(as.character(got), as.character(want)))
  cat(sprintf("  CHECK %-4s %-55s table %s | text %s\n", if (ok) "OK" else "DIFF", label,
              paste(got, collapse = " "), paste(want, collapse = " ")))
}
pct <- function(k, n, d = 1) formatC(100 * k / n, format = "f", digits = d)
sp  <- function(x) gsub("_", " ", x)

out  <- list()  # table name -> data.table, written at the end
NOTE <- list()  # table name -> footnote, written to table_notes.tsv
TITLE <- c(TableS1 = "Table S1. Accessions of the genome assemblies and reads used",
           TableS2 = "Table S2. BUSCO completeness of the BRAKER3/GALBA and ANNEVO annotations",
           TableS3 = "Table S3. MCMCtree priors and calibrations, with the effective prior on crown Anisoptera",
           TableS4 = "Table S4. Gene fate under three classifications (annotation only, TOGA2-validated, TOGA2 only)",
           TableS5 = "Table S5. Microchromosome orthogroup sharing: enrichment over chance and retention")
cat("repository root:", normalizePath(RUNS, mustWork = FALSE), "\n")

# Table S1: accessions of the assemblies and reads
cat("\n== Table S1: accessions ==\n")
if (have("assemblies", "kmer", "annot_sbatch")) {
  # one row per assembly, with what it was used for
  asm <- fread(p(IN$assemblies))                                    # genome, accession
  TEN   <- c(CODES, "Ischnura_elegans")
  MICRO <- c("Aeshna_juncea", "Ischnura_elegans", "Platycnemis_pennipes", "Tanypteryx_hageni")
  S1 <- asm[, .(data = "genome assembly", species = sp(genome), accession,
                details = "chromosome-level assembly", genome)]
  S1[genome == "Aeshna_juncea", details := "male, chromosome-level primary haplotype (ioAesJunc2.hap1.1)"]
  S1[genome == "Ischnura_elegans", details := "chromosome-level assembly; NCBI RefSeq annotation (ioIscEleg1.2) used for the ten-genome set"]
  for (i in seq_len(nrow(S1))) {
    g <- S1$genome[i]
    use <- "29-genome macrosynteny, ANNEVO annotation (Fig. 1, Table S2)"
    if (g %in% TEN)   use <- c(use, "ten-genome annotation and synteny (Fig. S1, Table S2)", "molecular dating (Fig. 3d)")
    if (g == "Aeshna_juncea") use <- c(use, "fusion junction, read depth, repeats, gene fate, k-mers (Figs 3, 4, S2-S4)")
    if (g == "Anax_junius")   use <- c(use, "junction synteny (Fig. 3c); gene-fate reference and TOGA2 reference (Fig. 4)")
    if (g == "Somatochlora_arctica") use <- c(use, "gene-fate outgroup (Fig. 4)")
    if (g %in% MICRO) use <- c(use, "microchromosome homology (Table S5)")
    S1[i, used_for := paste(use, collapse = "; ")]
  }

  # Ae. juncea hap2 and female assemblies and the HiFi reads of both individuals, as kmer.sh uses them
  # (kmer.sh names the hap2 fasta but not its accession)
  HAP2 <- "GCA_977065665.1"
  k <- readLines(p(IN$kmer), warn = FALSE)
  grab <- function(re) { x <- sub(re, "\\1", grep(re, k, value = TRUE, perl = TRUE), perl = TRUE); unique(x)[1] }
  fem_asm   <- grab("^ACC_FEM=(GCA_[0-9]+\\.[0-9]+).*")
  male_hifi <- grab("^MALE_READS=.*/([A-Z]{3}[0-9]+)\\.fastq.*")
  fem_hifi  <- grab("^FEM_READS=.*/([A-Z]{3}[0-9]+)\\.fastq.*")
  S1 <- rbind(S1, data.table(
    data      = c("genome assembly", "genome assembly", "HiFi reads", "HiFi reads"),
    species   = "Aeshna juncea",
    accession = c(HAP2, fem_asm, male_hifi, fem_hifi),
    details   = c("male, contig-level alternate haplotype (ioAesJunc2.hap2)",
                  "female, contig-level (isolate AK7, Alaska)",
                  "male (Scotland), the individual assembled",
                  "female (isolate AK7, Alaska; Tolman et al. 2025)"),
    used_for  = c("k-mer partition and alignment of the microchromosome and neo-Y (Supplementary Methods)",
                  "reciprocal k-mer analysis and alignment to the male assembly (Supplementary Methods)",
                  "read depth (Figs 3a, S2), neo-Y depth check, k-mers (Fig. S3) and heterozygosity",
                  "read depth (Figs 3a, S2), k-mers (Fig. S3) and heterozygosity"),
    genome = "Aeshna_juncea"), fill = TRUE)

  # RNA-seq used for annotation: the ACC<n>="..." assignments in each annotation script
  rna <- list()
  for (cd in names(CODES)) {
    x <- readLines(p(sprintf("annotation/genes/%s/%s.gene.annot.sbatch", CODES[[cd]], cd)), warn = FALSE)
    acc <- unique(sub('^ACC[0-9]*="([A-Z]{3}[0-9]+)".*', "\\1", grep('^ACC[0-9]*="[A-Z]{3}[0-9]+"', x, value = TRUE)))
    rna[[cd]] <- data.table(data = "RNA-seq reads", species = sp(CODES[[cd]]), accession = acc,
                            details = "RNA evidence for BRAKER3", used_for = "gene annotation (Fig. S1, Table S2)",
                            genome = CODES[[cd]])
  }
  rna <- rbindlist(rna)
  S1 <- rbind(S1, rna)

  # order: assemblies, reads, RNA-seq; by species, keeping the order above within a species
  S1[, `:=`(ord = match(data, c("genome assembly", "HiFi reads", "RNA-seq reads")), k = .I)]
  setorder(S1, ord, species, k)                                   # Ae. juncea: hap1, hap2, female
  S1[, c("ord", "k", "genome") := NULL]
  setnames(S1, c("Data", "Species", "Accession", "Details", "Used for"))
  print(S1[, .N, by = Data])
  CHECK("assemblies in the 29-genome set", nrow(asm), 29)
  CHECK("Ae. juncea male reads / female reads", c(male_hifi, fem_hifi), c("ERR12340099", "SRR36191036"))
  CHECK("species with RNA-seq", uniqueN(rna$species), 9)
  NOTE$TableS1 <- "Genome assemblies are the first (primary) haplotype of each species unless stated."
  out$TableS1 <- S1
}

# Table S2: BUSCO completeness, ANNEVO (29 genomes) and BRAKER3/GALBA (nine de-novo genomes)
cat("\n== Table S2: BUSCO ==\n")
if (have("assemblies", "annevo_busco", "annevo_genes", "braker_busco")) {
  asm <- fread(p(IN$assemblies))
  EPH <- c("Baetis", "Cloeon", "Ecdyonurus", "Ephemera", "Parafronurus", "Procloeon", "Serratella", "Siphlonurus")

  # ANNEVO BUSCO, one row per haplotype: keep the first haplotypes (_1.1 or no suffix; the _2.1 rows are not in
  # Fig. 1) and fix the one name misspelt in the csv
  ab <- fread(p(IN$annevo_busco))
  ab <- ab[!grepl("_2\\.1$", species)]
  ab[, genome := sub("_1\\.1$", "", species)]
  ab[genome == "Siphonurus_alternatus", genome := "Siphlonurus_alternatus"]

  # ANNEVO gene counts (genome with dots, genes, ...), same haplotype rule
  gc <- fread(p(IN$annevo_genes), header = FALSE)
  gc <- gc[!grepl("\\.2\\.1$", V1), .(genome = gsub(".", "_", V1, fixed = TRUE), genes = V2)]

  # BRAKER3/GALBA BUSCO short summaries: the count on the line holding each label ("<n>  Complete BUSCOs (C)" etc.)
  busco_n <- function(x, lab) as.integer(sub("^\\s*([0-9]+)\\s.*", "\\1", grep(lab, x, value = TRUE, fixed = TRUE)[1]))
  bb <- list()
  for (cd in names(CODES)) {
    x <- readLines(p(sprintf("annotation/genes/%s/busco_combined/%s", CODES[[cd]], BUSCO_TXT)), warn = FALSE)
    bb[[cd]] <- data.table(genome = CODES[[cd]], bC = busco_n(x, "(C)"), bS = busco_n(x, "(S)"), bD = busco_n(x, "(D)"),
                           bF = busco_n(x, "(F)"), bM = busco_n(x, "(M)"), bn = busco_n(x, "Total BUSCO groups searched"))
  }
  bb <- rbindlist(bb)

  # one row per assembly, Odonata first
  S2 <- merge(asm[, .(genome, accession)], ab[, .(genome, C, S, D, F, M, n)], by = "genome", all.x = TRUE)
  S2 <- merge(S2, gc, by = "genome", all.x = TRUE)
  S2 <- merge(S2, bb, by = "genome", all.x = TRUE)
  stopifnot(!anyNA(S2$C), !anyNA(S2$genes), nrow(S2) == nrow(asm))
  S2[, order := fifelse(sub("_.*", "", genome) %in% EPH, "Ephemeroptera", "Odonata")]
  setorder(S2, -order, genome)
  S2 <- S2[, .(Order = order, Species = sp(genome), Assembly = accession, `ANNEVO genes` = genes,
               `ANNEVO C (%)` = pct(C, n), `ANNEVO S (%)` = pct(S, n), `ANNEVO D (%)` = pct(D, n),
               `ANNEVO F (%)` = pct(F, n), `ANNEVO M (%)` = pct(M, n),
               `BRAKER3/GALBA C (%)` = fifelse(is.na(bn), NA_character_, pct(bC, bn)),
               `BRAKER3/GALBA S (%)` = fifelse(is.na(bn), NA_character_, pct(bS, bn)),
               `BRAKER3/GALBA D (%)` = fifelse(is.na(bn), NA_character_, pct(bD, bn)),
               `BRAKER3/GALBA F (%)` = fifelse(is.na(bn), NA_character_, pct(bF, bn)),
               `BRAKER3/GALBA M (%)` = fifelse(is.na(bn), NA_character_, pct(bM, bn)))]
  print(S2[, .(Species, `ANNEVO C (%)`, `BRAKER3/GALBA C (%)`)], nrows = 40)
  m <- ab[, .(C = mean(C / n)), by = .(o = fifelse(sub("_.*", "", genome) %in% EPH, "E", "O"))]
  CHECK("genomes per order (Odonata, Ephemeroptera)", S2[, .N, by = Order][order(-N), N], c(19, 10))
  CHECK("mean ANNEVO BUSCO C, Odonata (%)",       pct(m[o == "O", C], 1), "95.6")
  CHECK("mean ANNEVO BUSCO C, Ephemeroptera (%)", pct(m[o == "E", C], 1), "95.4")
  CHECK("ANNEVO genes per assembly (range)", format(range(gc[genome %in% asm$genome, genes]), big.mark = ","),
        c("13,720", "31,311"))
  CHECK("BUSCO groups (n), all runs", unique(c(ab$n, bb$bn)), 2708)
  NOTE$TableS2 <- paste("BUSCO in protein mode against insecta_odb12.2 (n = 2,708): C complete, S single-copy,",
                        "D duplicated, F fragmented, M missing. ANNEVO for all 29 genomes; BRAKER3/GALBA for the nine",
                        "de-novo annotations (Ischnura elegans uses its NCBI RefSeq annotation).")
  out$TableS2 <- S2
}

# Table S3: calibrations and effective priors, one row per node; the run settings and the priors on the
# birth-death process and the rates go in the table note
cat("\n== Table S3: MCMCtree ==\n")
if (have("mcmc_ctl", "calib_tree", "ages")) {
  # settings from the control file of the full-data run ("key = value * comment" lines)
  ctl <- trimws(sub("\\*.*$", "", readLines(p(IN$mcmc_ctl), warn = FALSE)))
  ctl <- ctl[grepl("=", ctl)]
  kv  <- setNames(trimws(sub("^[^=]*=", "", ctl)), trimws(sub("=.*$", "", ctl)))
  stopifnot(all(c("clock", "usedata", "model", "ncatG", "BDparas", "rgene_gamma", "sigma2_gamma",
                  "burnin", "sampfreq", "nsample") %in% names(kv)))
  num <- function(k) { v <- strsplit(kv[[k]], "[ \t]+")[[1]]; as.numeric(v[grepl("^[0-9.eE+-]+$", v)]) }   # drops the "C" of BDparas
  bd <- num("BDparas"); rg <- num("rgene_gamma"); s2 <- num("sigma2_gamma")
  gens <- num("burnin") + num("nsample") * num("sampfreq")
  cat(sprintf("  ctl: clock %s, usedata %s, model %s, ncatG %s, BDparas %s, rgene_gamma %s, sigma2_gamma %s\n",
              kv[["clock"]], kv[["usedata"]], kv[["model"]], kv[["ncatG"]], kv[["BDparas"]], kv[["rgene_gamma"]], kv[["sigma2_gamma"]]))
  stopifnot(kv[["clock"]] == "2", kv[["model"]] == "7", grepl("^2", kv[["usedata"]]))   # what the note says
  CHECK("generations (burn-in + samples x frequency)", gens, 2.2e6)

  # calibrated tree: each 'B(min, max, pL, pU)' label becomes node label CAL<i> so ape can read it; ape numbers
  # the internal nodes as MCMCtree's t_n (root = tips + 1), which the CHECKs below confirm
  tre <- paste(readLines(p(IN$calib_tree), warn = FALSE)[-1], collapse = "")      # first line is "11 1"
  bnd <- regmatches(tre, gregexpr("'B\\([^)]*\\)'", tre))[[1]]
  for (i in seq_along(bnd)) tre <- sub(bnd[i], paste0("CAL", i), tre, fixed = TRUE)
  tr <- read.tree(text = tre)
  nt <- Ntip(tr)

  # clade names and fossils (Methods), keyed by the tips each clade contains
  ANI <- c("Aeshna_juncea_neoX", "Aeshna_juncea_neoY", "Anax_junius", "Cordulegaster_boltonii",
           "Somatochlora_arctica", "Pantala_flavescens", "Orthetrum_coerulescens", "Sympetrum_striolatum")
  CLADES <- list(
    list("crown Odonata",                  "Triassolestodes",   tr$tip.label),
    list("crown Anisoptera",               "Sinacymatophlebia", ANI),
    list("Cavilabiata",                    "Juralibellula",     ANI[4:8]),
    list("Libellulidae + Corduliidae",     "Palaeolibellula",   ANI[5:8]),
    list("crown Libellulidae",             "Tauriphila",        ANI[6:8]),
    list("crown Zygoptera",                "Mersituria",        c("Ischnura_elegans", "Platycnemis_pennipes", "Archilestes_grandis")),
    list("Aeshna + Anax",                  NA,                  ANI[1:3]),
    list("neo-X / neo-Y gametologs",       NA,                  ANI[1:2]),
    list("Pantala + Orthetrum",            NA,                  ANI[6:7]),
    list("Ischnura + Platycnemis",         NA,                  c("Ischnura_elegans", "Platycnemis_pennipes")))
  key <- function(x) paste(sort(x), collapse = ",")
  ck  <- sapply(CLADES, function(cl) key(cl[[3]]))

  # one row per internal node: its clade (or tip list), fossil and calibration bounds
  S3 <- list()
  for (j in seq_len(tr$Nnode)) {
    nd   <- nt + j
    tips <- extract.clade(tr, nd)$tip.label
    i    <- match(key(tips), ck)
    lab  <- tr$node.label[j]
    b    <- if (grepl("^CAL", lab)) gsub("'", "", bnd[as.integer(sub("CAL", "", lab))]) else NA_character_
    bv   <- if (is.na(b)) rep(NA_real_, 4) else as.numeric(strsplit(gsub("[B() ]", "", b), ",")[[1]])
    S3[[j]] <- data.table(node = paste0("t_n", nd), clade = if (is.na(i)) key(tips) else CLADES[[i]][[1]],
                          fossil = if (is.na(i)) NA_character_ else CLADES[[i]][[2]], calibration = b,
                          min_Mya = 100 * bv[1], max_Mya = 100 * bv[2])
  }
  S3 <- rbindlist(S3)

  # effective priors: the prior-only run (time unit 100 Myr)
  ag <- fread(p(IN$ages))[run == "prior_only"]
  S3 <- merge(S3, ag[, .(node, mean, et_lo, et_hi, hpd_lo, hpd_hi)], by = "node", all.x = TRUE, sort = FALSE)
  stopifnot(!anyNA(S3$mean))
  f1 <- function(x) formatC(100 * x, format = "f", digits = 1)
  S3 <- S3[, .(Node = node, Clade = clade, Fossil = fossil, Calibration = calibration,
               `Minimum (Mya)` = min_Mya, `Maximum (Mya)` = max_Mya,
               `Prior mean (Mya)` = f1(mean), `Prior 2.5-97.5% (Mya)` = paste0(f1(et_lo), "-", f1(et_hi)),
               `Prior 95% HPD (Mya)` = paste0(f1(hpd_lo), "-", f1(hpd_hi)))]
  print(S3)
  CHECK("node of the neo-X/neo-Y gametologs", S3[Clade == "neo-X / neo-Y gametologs", Node], "t_n15")
  CHECK("calibrated nodes", S3[!is.na(Calibration), .N], 6)
  CHECK("crown Anisoptera, effective prior 2.5% (Mya)",
        round(100 * ag[node == S3[Clade == "crown Anisoptera", Node], et_lo]), 181)
  CHECK("gametolog node, prior mean and 95% HPD (Mya)",
        S3[Node == "t_n15", c(`Prior mean (Mya)`, `Prior 95% HPD (Mya)`)], c("85.3", "0.0-223.4"))
  NOTE$TableS3 <- sprintf(paste(
    "Ages in Mya (MCMCtree time unit 100 Myr). Calibrations are soft-bounded uniform priors B(min, max, pL, pU)",
    "with 2.5%% prior probability beyond each bound; effective priors are the node ages sampled without data",
    "(prior-only run). All runs used an independent-rates relaxed clock with the approximate likelihood",
    "(Hessian from each dataset, GTR with %s gamma categories), a birth-death prior on node ages",
    "(birth rate %g, death rate %g, sampling fraction %g), gamma-Dirichlet priors on the mean rate",
    "(G(%g, %g), mean %g substitutions per site per 100 Myr) and on the rate variance among branches",
    "(G(%g, %g)), and %s iterations (%s burn-in, then %s samples every %s)."),
    kv[["ncatG"]], bd[1], bd[2], bd[3], rg[1], rg[2], rg[1] / rg[2], s2[1], s2[2],
    format(gens, big.mark = ",", scientific = FALSE), format(num("burnin"), big.mark = ",", scientific = FALSE),
    format(num("nsample"), big.mark = ",", scientific = FALSE), kv[["sampfreq"]])
  out$TableS3 <- S3
}

# Table S4: gene fate under three classifications, from figs/fig4/fig4.R; without fig4.R, a TableS4.tsv
# in this folder is read back, checked and written out with the rest
cat("\n== Table S4: gene fate ==\n")
S4 <- NULL
if (file.exists(FIG4)) {
  # run fig4.R in its own folder and environment, output hidden (it rewrites its own outputs there)
  F4 <- new.env()
  invisible(capture.output(ok4 <- tryCatch({ sys.source(FIG4, envir = F4, chdir = TRUE); TRUE },
    error = function(e) { message("  fig4.R failed (run it alone to see where): ", conditionMessage(e)); FALSE })))
  if (ok4) {
    CLS  <- c("autosomes", "ancestral X", "neo-X", "neo-Y")
    CATS <- setdiff(F4$segs, "excluded")

    # the classifications: count matrix (class x category), categories and their descriptions
    SETS <- list(`annotation only` = list(M = F4$fig_sets$annot$M, cats = CATS, desc = F4$ctitle[CATS]),
                 `TOGA2-validated` = list(M = F4$fig_sets$toga2$M, cats = CATS, desc = F4$ctitle[CATS]))
    if (!is.null(F4$tt)) {
      SETS$`TOGA2 only` <- list(M = F4$tt, cats = colnames(F4$tt),
                                desc = c("intact projection on its own chromosome class", "intact only on the sister arm (= degenerated)",
                                         "intact only elsewhere in the genome", "no intact projection"))
    }

    # per classification: n and "% (95% CI)" (Wilson, as in fig4.R) per class, then a row of denominators
    for (label in names(SETS)) {
      M <- SETS[[label]]$M[CLS, SETS[[label]]$cats, drop = FALSE]; tot <- rowSums(M)
      b <- data.table(Classification = label, Category = SETS[[label]]$cats, Description = unname(SETS[[label]]$desc))
      for (cl in CLS) {
        txt <- character()
        for (k in M[cl, ]) {
          ci  <- F4$wilson(k, tot[[cl]])
          txt <- c(txt, sprintf("%s (%s-%s)", pct(k, tot[[cl]]), pct(ci[1], 1), pct(ci[2], 1)))
        }
        b[, (paste(cl, "n")) := as.integer(M[cl, ])]
        b[, (paste(cl, "% (95% CI)")) := txt]
      }
      inf <- data.table(Classification = label, Category = "informative genes", Description = "denominator")
      for (cl in CLS) inf[, c(paste(cl, "n"), paste(cl, "% (95% CI)")) := .(as.integer(tot[[cl]]), "")]
      S4 <- rbind(S4, b, inf)
    }
    if (is.null(F4$tt)) message("  fig4.R found no TOGA2 projections: TOGA2-only block left out")
  }
} else {
  f4 <- c("TableS4.tsv", "TableS4a.tsv")[file.exists(c("TableS4.tsv", "TableS4a.tsv"))][1]   # S4a: earlier name
  if (is.na(f4)) message("  no ../figs/fig4/fig4.R and no TableS4.tsv (skipped)")
  else { S4 <- fread(f4); cat("  read", f4, "\n") }
}
# print the neo-Y columns and check them against the text; g() = a neo-Y cell of the table
if (!is.null(S4)) {
  print(S4[, c(1, 2, grep("neo-Y", names(S4))), with = FALSE])
  g <- function(cl, cat, col) S4[Classification == cl & Category == cat, get(paste("neo-Y", col))]
  CHECK("neo-Y conserved, annotation only (%)",        sub(" .*", "", g("annotation only", "conserved", "% (95% CI)")), "24.1")
  CHECK("neo-Y degenerated, annotation only (95% CI)", g("annotation only", "degenerated", "% (95% CI)"), "65.5 (60.9-69.8)")
  CHECK("neo-Y conserved / degenerated, TOGA2-validated (%)",
        sub(" .*", "", c(g("TOGA2-validated", "conserved", "% (95% CI)"), g("TOGA2-validated", "degenerated", "% (95% CI)"))),
        c("55.6", "40.3"))
  CHECK("neo-X conserved / degenerated, annotation only (%)",
        sub(" .*", "", S4[Classification == "annotation only" & Category %in% c("conserved", "degenerated"), `neo-X % (95% CI)`]),
        c("86.3", "0.9"))
  CHECK("neo-Y informative genes, TOGA2 only", g("TOGA2 only", "informative genes", "n"), 432)
  CHECK("neo-Y intact on arm / sister only, TOGA2 only (%)",
        sub(" .*", "", c(g("TOGA2 only", "intact on arm", "% (95% CI)"), g("TOGA2 only", "intact on sister only", "% (95% CI)"))),
        c("51.2", "42.6"))
  NOTE$TableS4 <- paste("Informative focal An. junius genes per chromosome class; % with Wilson 95% CI. TOGA2 only",
                        "scores presence of an intact, orthologous, non-retro projection, not copy number.")
  out$TableS4 <- S4
}

# Table S5: microchromosome orthogroup sharing, one row per species pair: shared orthogroups against chance,
# and retention in each direction (same counts as microchromosome/ortholog_sharing/mchrom_pairwise.py)
cat("\n== Table S5: microchromosome sharing ==\n")
if (have("n0", "micro_beds", "juncea_cds")) {
  MICRO_CHR <- c(Aeshna_juncea = "OZ369553.1", Ischnura_elegans = "NC_060258.1",
                 Platycnemis_pennipes = "OW121858.1", Tanypteryx_hageni = "CM054119.1")
  BEDS <- c(Aeshna_juncea = IN$juncea_cds, setNames(IN$micro_beds, names(MICRO_CHR)[-1]))
  norm <- function(x) sub("\\.t[0-9]+$", "", trimws(x))           # the --strip-re '\.t\d+$' of the run
  n0 <- fread(p(IN$n0), sep = "\t", select = c("HOG", names(MICRO_CHR)), colClasses = "character")
  total <- nrow(n0)

  # per species: the HOGs with a gene on its microchromosome, and the HOGs it has anywhere
  on_m <- list(); any_p <- list()
  for (s in names(MICRO_CHR)) {
    bed <- fread(p(BEDS[[s]]), header = FALSE, select = c(1, 4), colClasses = "character")
    onm_genes <- unique(norm(bed[V1 == MICRO_CHR[[s]], V4]))
    stopifnot(length(onm_genes) > 0)
    z <- n0[nzchar(trimws(get(s))), .(HOG, cell = get(s))]
    x <- z[, .(gene = norm(unlist(strsplit(cell, ",", fixed = TRUE)))), by = HOG]
    any_p[[s]] <- unique(x$HOG)
    on_m[[s]]  <- unique(x[gene %chin% onm_genes, HOG])
  }

  # retention a -> b: of the HOGs on a's microchromosome that b carries anywhere, the share on b's
  ret <- function(a, b) {
    carried <- intersect(on_m[[a]], any_p[[b]]); k <- length(intersect(carried, on_m[[b]]))
    list(k = k, n = length(carried), txt = sprintf("%s (%d/%d)", pct(k, length(carried)), k, length(carried)))
  }

  # one row per species pair; rr collects the retention fractions, both directions
  S5 <- list(); rr <- c()
  for (pr in combn(names(MICRO_CHR), 2, simplify = FALSE)) {
    a <- pr[1]; b <- pr[2]
    obs <- length(intersect(on_m[[a]], on_m[[b]])); ex <- length(on_m[[a]]) * length(on_m[[b]]) / total
    S5[[paste(a, b)]] <- data.table(`Species A` = sp(a), `Species B` = sp(b),
               `HOGs on A's microchromosome` = length(on_m[[a]]), `HOGs on B's microchromosome` = length(on_m[[b]]),
               `Shared (observed)` = obs, `Expected if independent` = round(ex, 1), `Fold enrichment` = round(obs / ex, 1),
               `Retention A to B, % (k/n)` = ret(a, b)$txt, `Retention B to A, % (k/n)` = ret(b, a)$txt)
    rr <- c(rr, ret(a, b)$k / ret(a, b)$n, ret(b, a)$k / ret(b, a)$n)
  }
  S5 <- rbindlist(S5)
  print(S5)
  CHECK("total HOGs", format(total, big.mark = ","), "32,481")
  CHECK("HOGs on each microchromosome (Aj, Ie, Pp, Th)", sapply(on_m, length), c(261, 246, 257, 133))
  CHECK("fold enrichment range", format(range(S5$`Fold enrichment`), nsmall = 1), c("10.6", "50.5"))
  CHECK("Ae. juncea - T. hageni shared / expected",
        unlist(S5[`Species A` == "Aeshna juncea" & `Species B` == "Tanypteryx hageni", .(`Shared (observed)`, `Expected if independent`)]),
        c(54, 1.1))
  CHECK("retention range (%)", as.character(round(100 * range(rr))), c("20", "66"))
  NOTE$TableS5 <- paste("Hierarchical orthogroups (HOGs, N0) of the 29-genome OrthoFinder run. Expected = HOGs on A's",
                        "microchromosome x HOGs on B's / all HOGs. Retention A to B: of the HOGs on A's microchromosome",
                        "that B carries anywhere (n), those on B's microchromosome (k). Microchromosomes:",
                        paste(sp(names(MICRO_CHR)), MICRO_CHR, collapse = "; "))
  out$TableS5 <- S5
}

# write one file per table, plus the titles and notes for the legends
xl <- requireNamespace("writexl", quietly = TRUE)
note <- character()
for (nm in names(out)) {
  fwrite(out[[nm]], paste0(nm, ".tsv"), sep = "\t", na = "")
  if (xl) writexl::write_xlsx(list(as.data.frame(out[[nm]])), paste0(nm, ".xlsx"))
  note <- c(note, if (is.null(NOTE[[nm]])) "" else NOTE[[nm]])
}
notes <- data.table(Table = names(out), Title = TITLE[names(out)], Note = note)
fwrite(notes, "table_notes.tsv", sep = "\t")
cat("\nwritten:", paste0(names(out), ".tsv"), if (xl) paste0(names(out), ".xlsx"), "table_notes.tsv\n")
if (!xl && length(out)) cat("(no writexl in this R, so TSV only; install.packages(\"writexl\") for .xlsx)\n")
