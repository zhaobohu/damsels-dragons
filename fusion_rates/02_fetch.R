#!/usr/bin/env Rscript
# Fetch from NCBI what BOLD lacks (COI for gap species, all other markers), batched and cached,
# then rebuild unfiltered.COI.fna, prep.csv and raw_<marker>.fasta from the cache

suppressPackageStartupMessages({ library(data.table); library(rentrez); library(Biostrings) })
source("config.R")

if (nzchar(ENTREZ_KEY)) set_entrez_key(ENTREZ_KEY)

if (!nzchar(ENTREZ_EMAIL)) stop("export ENTREZ_EMAIL=<your address> (read in config.R)", call. = FALSE)

# pause between calls (NCBI rate limit: 10/s with a key, 3/s without)
pause <- function() Sys.sleep(if (nzchar(ENTREZ_KEY)) 0.11 else 0.34)

# length inside the marker's MARKER_LEN window (MIN_LEN-MAX_LEN if not listed)
len_ok <- function(marker, len) {
  bnd <- function(m, i) vapply(m, function(x) {
    r <- MARKER_LEN[[x]]; if (is.null(r)) c(MIN_LEN, MAX_LEN)[i] else r[i] }, 0)
  !is.na(len) & len >= bnd(marker, 1L) & len <= bnd(marker, 2L)
}

# gene-name query terms per marker
alias <- list(
  COI  = '(COI[All Fields] OR COX1[All Fields] OR "cytochrome c oxidase subunit I"[All Fields])',
  EF1a = '(EF-1a[All Fields] OR EF1alpha[All Fields] OR "elongation factor 1"[All Fields])')

# inputs; BOLD accessions are skipped as duplicates
taxa <- fread(TAXA_TSV)
prep <- fread("prep.csv")
gap  <- readLines("needed_from_ncbi.txt"); gap <- gap[nzchar(gap)]
seen <- unique(prep$insdc_acs[nzchar(as.character(prep$insdc_acs))])

# cache of every sequence downloaded so far; records outside the current windows are re-queried
cache <- if (file.exists("ncbi_cache.csv")) {
  fread("ncbi_cache.csv", colClasses = list(character = "seq"))
} else {
  data.table(species = character(), marker = character(),
             accession = character(), seq = character(), len = integer())
}
if (nrow(cache)) {
  cache[, len := nchar(seq)]
  n0 <- nrow(cache)
  cache <- cache[len_ok(marker, len)]
  if (nrow(cache) < n0)
    cat(sprintf("cache: dropped %d records failing the current MARKER_LEN windows\n",
                n0 - nrow(cache)))
}

# cache contents, for this run's markers and for others (ignored)
cat(sprintf("markers in this run: %s\n", paste(MARKERS, collapse = ", ")))
cat(sprintf("cache: %d sequences total, %d for these markers (%d species)\n",
            nrow(cache), nrow(cache[marker %chin% MARKERS]),
            uniqueN(cache[marker %chin% MARKERS]$species)))
if (nrow(cache[!marker %chin% MARKERS]))
  cat(sprintf("       %d cached sequences are for other markers (%s) and are ignored\n",
              nrow(cache[!marker %chin% MARKERS]),
              paste(sort(unique(cache[!marker %chin% MARKERS]$marker)), collapse = ", ")))

# species-marker pairs still to download (COI only for species without BOLD records)
want <- rbindlist(lapply(MARKERS, function(m)
  data.table(species = if (m == "COI") gap else taxa$species, marker = m)))
want <- unique(want[nzchar(species)])
want[, n := fifelse(marker == "COI", N_COI, N_MARKER)]
need <- want[!cache[, .(species, marker)], on = .(species, marker)]
cat(sprintf("to download: %d species-marker pairs (%d already cached)\n",
            nrow(need), nrow(want) - nrow(need)))

# binomial from an NCBI organism name (drops subspecies and strain suffixes)
first2 <- function(x) {
  p <- strsplit(x, " ", fixed = TRUE)
  vapply(p, function(v) if (length(v) >= 2) paste(v[1], v[2]) else NA_character_, "")
}

# say why a fetch returned nothing
.why <- function(...) cat("      ", sprintf(...), "\n", sep = "")

# first line of an error (NCBI sends whole HTML pages)
.short <- function(e) substr(sub("\n.*", "", conditionMessage(e)), 1, 140)

# retry transient NCBI errors with exponential backoff (expr is re-evaluated each try)
.TRANSIENT <- paste0("HTTP failure:? *(429|500|502|503|504)|Timeout|",
                     "Could not resolve|Connection reset|Operation timed out")
retry <- function(expr, what, verbose = TRUE) {
  call <- substitute(expr); env <- parent.frame()
  for (k in seq_len(RETRY_TRIES)) {
    r <- tryCatch(eval(call, env), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    if (!grepl(.TRANSIENT, conditionMessage(r)) || k == RETRY_TRIES) return(r)
    wait <- RETRY_WAIT * 2^(k - 1)
    if (verbose) .why("%s failed (%s) -- retry %d/%d in %.0fs",
                      what, .short(r), k, RETRY_TRIES - 1, wait)
    Sys.sleep(wait)
  }
}

# one batch: esearch an OR-list of species, page esummary, map records to species, efetch;
# returns a fetch_error object on failure so the caller can bisect
fetch_batch <- function(spp, marker, n, verbose = TRUE) {
  rng <- MARKER_LEN[[marker]]; if (is.null(rng)) rng <- c(MIN_LEN, MAX_LEN)
  term <- sprintf("(%s) AND %s AND %d:%d[SLEN]",
                  paste(sprintf('"%s"[Organism]', spp), collapse = " OR "),
                  alias[[marker]], rng[1], rng[2])
  s <- retry(entrez_search("nuccore", term, retmax = 0, use_history = TRUE,
                             email = ENTREZ_EMAIL, tool = ENTREZ_TOOL),
             "esearch", verbose = verbose)
  pause()
  if (inherits(s, "error")) {
    if (verbose) .why("esearch FAILED (term %d chars): %s", nchar(term), .short(s))
    return(structure(list(), class = "fetch_error"))
  }
  if (s$count == 0) {
    if (verbose) .why("esearch returned 0 hits (term %d chars)", nchar(term))
    return(NULL)
  }
  if (is.null(s$web_history)) {
    if (verbose) .why("esearch gave %d hits but NO web_history -- use_history failed",
                      s$count)
    return(structure(list(), class = "fetch_error"))
  }

  # esummary, 300 at a time
  recs <- list()
  for (st in seq(0, s$count - 1, by = 300)) {
    su <- retry(entrez_summary("nuccore", web_history = s$web_history,
                                email = ENTREZ_EMAIL, tool = ENTREZ_TOOL,
                               retstart = st, retmax = 300),
                "esummary", verbose = verbose)
    pause()
    if (inherits(su, "error")) { if (verbose) .why("esummary gave up: %s", .short(su)); next }
    if (!is.null(su$uid)) su <- list(su)          # a single hit is not a list
    recs <- c(recs, su)
  }
  if (!length(recs)) {
    if (verbose) .why("esearch found %d but esummary returned nothing", s$count)
    return(structure(list(), class = "fetch_error"))
  }

  # uid, accession, organism and length of each record
  g <- function(x, f) { v <- x[[f]]; if (is.null(v)) NA else v }
  d <- data.table(
    uid = vapply(recs, function(x) as.character(g(x, "uid")), ""),
    acc = vapply(recs, function(x) as.character(g(x, "caption")), ""),
    org = vapply(recs, function(x) as.character(g(x, "organism")), ""),
    len = vapply(recs, function(x) as.integer(g(x, "slen")), 0L))

  # NCBI organism name, collapsed to the binomial
  if (all(is.na(d$org)) && verbose)
    .why("esummary has no 'organism' field -- fields seen: %s",
         paste(head(names(recs[[1]]), 20), collapse = ", "))
  d[, ncbi_name := first2(org)]
  d0org <- d$ncbi_name
  n0 <- nrow(d)

  # one-species query: everything belongs to that species (NCBI may use a synonym), except
  # records NCBI assigns to another of our species; new NCBI names are learned as aliases
  if (length(spp) == 1L) {
    steal <- !is.na(d$ncbi_name) & d$ncbi_name != spp & d$ncbi_name %chin% .ALLSPP
    if (any(steal)) {
      .BLOCKED[[spp]] <<- unique(c(.BLOCKED[[spp]], unique(d$ncbi_name[steal])))
      if (verbose) .why("dropped %d records NCBI attributes to %s, not %s",
                        sum(steal), paste(unique(d$ncbi_name[steal]), collapse = ", "), spp)
      d <- d[!steal]
      if (!nrow(d)) return(NULL)
    }
    d[, species := spp]
    newal <- unique(d[!is.na(ncbi_name) & ncbi_name != spp, ncbi_name])
    for (a in newal) {
      prev <- .ALIAS[[a]]
      if (!is.null(prev) && prev != spp) {
        # two of our species share this NCBI name: retire the alias
        .AMBIG[[a]] <<- TRUE
        .ALIAS[[a]] <<- NULL
        if (verbose) .why("alias '%s' is AMBIGUOUS (%s and %s) -- not used in batches",
                          a, prev, spp)
      } else if (is.null(.AMBIG[[a]])) {
        .ALIAS[[a]] <<- spp
        if (verbose) .why("alias: NCBI calls %s '%s'", spp, a)
      }
    }
  } else {
    # batch: match on the organism name, then on learned aliases
    d[, species := ncbi_name]
    al  <- unlist(as.list(.ALIAS))
    hit <- !d$species %chin% spp & d$ncbi_name %chin% names(al)
    if (any(hit)) d[hit, species := unname(al[ncbi_name])]
  }

  # requested species only, no BOLD duplicates, inside the length window; n longest per species
  d <- d[species %chin% spp]
  n1 <- nrow(d)
  d <- d[!acc %chin% seen & len_ok(marker, len)]
  if (!nrow(d)) {
    if (verbose) {
      .why("%d summaries -> %d matched a requested species -> 0 after dedupe/length", n0, n1)
      if (n0 > 0 && n1 == 0 && length(spp) > 1)
        .why("  NCBI returned these organisms: %s",
             paste(head(unique(na.omit(d0org)), 6), collapse = "; "))
    }
    return(NULL)
  }
  setorder(d, species, -len)
  d <- d[, head(.SD, n), by = species]

  # efetch, 200 ids at a time
  seqs <- list()
  for (i in seq(1, nrow(d), by = 200)) {
    ids <- d$uid[i:min(i + 199, nrow(d))]
    fa <- retry(entrez_fetch("nuccore", id = ids, rettype = "fasta",
                          email = ENTREZ_EMAIL, tool = ENTREZ_TOOL),
                "efetch", verbose = verbose)
    pause()
    if (inherits(fa, "error")) { if (verbose) .why("efetch gave up: %s", .short(fa)); next }
    blocks <- strsplit(sub("^>", "", fa), "\n>", fixed = TRUE)[[1]]
    for (b in blocks) {
      ln <- strsplit(b, "\n", fixed = TRUE)[[1]]
      if (length(ln) < 2) next
      a <- sub("\\..*$", "", sub("^(\\S+).*$", "\\1", ln[1]))   # accession, no version
      sq <- paste(ln[-1], collapse = "")
      if (nchar(sq)) seqs[[a]] <- sq
    }
  }
  # sequence per accession, NA where efetch returned none (unlist(seqs[acc]) would shift rows)
  d[, seq := vapply(acc, function(a) {
    v <- seqs[[a]]; if (is.null(v)) NA_character_ else v }, "")]
  d <- d[!is.na(seq) & nchar(seq) > 0]
  if (!nrow(d)) return(NULL)
  d[, .(species, marker = marker, accession = acc, seq, len = nchar(seq))]
}

# split species into chunks whose esearch term stays under MAX_TERM_CHARS (at most BATCH_SIZE each)
chunk_by_term <- function(spp, marker) {
  rng <- MARKER_LEN[[marker]]; if (is.null(rng)) rng <- c(MIN_LEN, MAX_LEN)
  fixed <- nchar(sprintf("() AND %s AND %d:%d[SLEN]", alias[[marker]], rng[1], rng[2]))
  cost  <- nchar(spp) + nchar('""[Organism] OR ')
  out <- list(); cur <- integer(0); tot <- fixed
  for (i in seq_along(spp)) {
    if (length(cur) && (tot + cost[i] > MAX_TERM_CHARS || length(cur) >= BATCH_SIZE)) {
      out[[length(out) + 1]] <- spp[cur]; cur <- integer(0); tot <- fixed
    }
    cur <- c(cur, i); tot <- tot + cost[i]
  }
  if (length(cur)) out[[length(out) + 1]] <- spp[cur]
  out
}

# NCBI-name -> our-name aliases (environments, so single aliases can be retired)
.ALLSPP <- unique(c(taxa$species, gap))   # every name we could attribute to
.BLOCKED <- new.env(parent = emptyenv())  # species whose hits belong to someone else
.ALIAS <- new.env(parent = emptyenv())
.AMBIG <- new.env(parent = emptyenv())
if (file.exists("ncbi_alias.csv")) {
  a <- fread("ncbi_alias.csv")
  for (i in seq_len(nrow(a))) assign(a$ncbi_name[i], a$species[i], envir = .ALIAS)
  cat(sprintf("aliases: %d NCBI names mapped\n", nrow(a)))
}

# download each marker in term-length chunks
new_rows <- list()
for (m in MARKERS) {
  spp <- need[marker == m, species]
  if (!length(spp)) next
  n <- need[marker == m, n[1]]
  # fetch a chunk, bisecting on failure down to single species (depth >= 2 hides "0 hits")
  try_chunk <- function(sp, depth = 0L) {
    r <- tryCatch(fetch_batch(sp, m, n, verbose = depth < 2L),
                  error = function(e) structure(list(), class = "fetch_error"))
    if (!inherits(r, "fetch_error")) return(r)
    if (length(sp) == 1L) {
      cat(sprintf("      giving up on %s\n", sp)); return(NULL)
    }
    h <- ceiling(length(sp) / 2)
    cat(sprintf("      batch of %d failed -- bisecting\n", length(sp)))
    rbindlist(list(try_chunk(sp[1:h], depth + 1L),
                   try_chunk(sp[(h + 1):length(sp)], depth + 1L)), fill = TRUE)
  }

  chunks <- chunk_by_term(spp, m)
  got <- 0L; hits <- list()
  for (ci in seq_along(chunks)) {
    r <- try_chunk(chunks[[ci]])
    if (!is.null(r) && nrow(r)) {
      new_rows[[length(new_rows) + 1]] <- r; hits[[length(hits) + 1]] <- r
      got <- got + nrow(r)
    }
    cat(sprintf("%-5s chunk %d/%d (%d species) -> %d records\n",
                m, ci, length(chunks), length(chunks[[ci]]), if (is.null(r)) 0L else nrow(r)))
  }
  found <- if (length(hits)) unique(rbindlist(hits, fill = TRUE)$species) else character(0)

  # sweep: re-query the species the batches missed, one per query
  # (RECOVERED ~0 means the batching is faithful; large means lower MAX_TERM_CHARS)
  miss <- setdiff(spp, found)
  rec <- 0L; rec_sp <- 0L
  if (length(miss) && SWEEP_SIZE > 0) {
    sw <- split(miss, ceiling(seq_along(miss) / SWEEP_SIZE))
    for (si in seq_along(sw)) {
      r <- try_chunk(sw[[si]], depth = 2L)
      if (!is.null(r) && nrow(r)) {
        new_rows[[length(new_rows) + 1]] <- r
        rec <- rec + nrow(r); rec_sp <- rec_sp + uniqueN(r$species)
      }
    }
    cat(sprintf("%-5s sweep: re-queried %d species the batch missed -> RECOVERED %d records for %d species\n",
                m, length(miss), rec, rec_sp))
  }
  cat(sprintf("%-5s total %d records for %d species (%d from batch, %d from sweep)\n",
              m, got + rec, length(found) + rec_sp, length(found), rec_sp))
}

# save learned aliases; report blocked species and ambiguous aliases
al <- unlist(as.list(.ALIAS))
if (length(al)) {
  fwrite(data.table(ncbi_name = names(al), species = unname(al)), "ncbi_alias.csv")
  cat(sprintf("\naliases: %d NCBI names mapped to our names (ncbi_alias.csv)\n", length(al)))
}
bl <- ls(.BLOCKED)
if (length(bl)) {
  cat(sprintf("\nBLOCKED -- NCBI returned records under a different species name for\n%d of our species, so nothing was attributed to them. Check these by hand:\n", length(bl)))
  for (b in bl) cat(sprintf("  %-32s -> NCBI returned %s\n", b,
                            paste(.BLOCKED[[b]], collapse = ", ")))
  cat("A search for a name NCBI does not recognise gets split by automatic term\nmapping and can match an unrelated species. Verify the name, or add a\nGBIF_OVERRIDE / PRE_FIX entry in reconcile.R.\n")
}

if (length(ls(.AMBIG)))
  cat(sprintf("aliases: %d NCBI names were ambiguous and are not used in batches: %s\n",
              length(ls(.AMBIG)), paste(head(ls(.AMBIG), 5), collapse = "; ")))

# add the new downloads to the cache
if (length(new_rows)) {
  add <- rbindlist(new_rows, fill = TRUE)
  add <- add[, .(species, marker, accession, seq, len)]
  cache <- unique(rbind(cache, add, fill = TRUE), by = c("species", "marker", "accession"))
  fwrite(cache, "ncbi_cache.csv")
  cat(sprintf("\ncache: added %d newly downloaded sequences (%d total)\n",
              nrow(add), nrow(cache)))
} else if (nrow(cache)) {
  fwrite(cache, "ncbi_cache.csv")
  cat("\ncache: nothing new downloaded\n")
}

# rebuild (never append) the COI pool: BOLD records from 01_prep.R plus cached NCBI COI
bold_ids <- prep[source == "bold", record_id]
allseq   <- Biostrings::readDNAStringSet("unfiltered.COI.fna")
names(allseq) <- sub("^_R_", "", names(allseq))
allseq   <- allseq[!duplicated(names(allseq))]
boldseq  <- allseq[names(allseq) %chin% bold_ids]

coi <- cache[marker == "COI" & species %chin% gap]
coi[, record_id := gsub("[^A-Za-z0-9]+", "_", accession)]
coi <- coi[!duplicated(record_id)]

out <- c(as.character(boldseq), setNames(coi$seq, coi$record_id))
writeLines(as.vector(rbind(sprintf(">%s", names(out)), unname(out))),
           "unfiltered.COI.fna")

# rebuild prep.csv the same way
add <- if (nrow(coi)) data.table(record_id = coi$record_id,
  gbif_key = taxa$gbif_key[match(coi$species, taxa$species)],
  acc_name = coi$species, species_bold = NA_character_,
  family_bold = taxa$family[match(coi$species, taxa$species)],
  bin_uri = NA_character_, insdc_acs = coi$accession, source = "ncbi") else NULL
newprep <- rbind(prep[source == "bold"], add, fill = TRUE)
stopifnot(!anyDuplicated(newprep$record_id))
fwrite(newprep, "prep.csv")
cat(sprintf("COI pool rebuilt: %d BOLD + %d NCBI records (%d NCBI species)\n",
            length(boldseq), nrow(coi), uniqueN(coi$species)))

# other markers: the longest cached record per species -> raw_<marker>.fasta
for (m in setdiff(MARKERS, "COI")) {
  d <- cache[marker == m & species %chin% taxa$species]
  if (!nrow(d)) next
  setorder(d, species, -len)
  d <- d[!duplicated(species)]
  writeLines(as.vector(rbind(sprintf(">%s", gsub(" ", "_", d$species)), d$seq)),
             sprintf("raw_%s.fasta", m))
  cat(sprintf("%-5s %d sequences, length %d-%d (median %d)\n", m, nrow(d),
              min(d$len), max(d$len), as.integer(median(d$len))))
  if (max(d$len) > 2 * median(d$len))
    cat(sprintf("  WARNING: %s has records >2x the median -- likely multi-gene\n", m))
}

# markers per species
per <- cache[species %chin% taxa$species, uniqueN(marker), by = species]
cat(sprintf("\n0 markers: %d | 1 marker: %d | >=3 markers: %d (of %d species)\n",
            nrow(taxa) - nrow(per), per[V1 == 1, .N], per[V1 >= 3, .N], nrow(taxa)))
