# Turns the stochastic maps into fusion events for 13_figures.R:
# n_iter_available() and events_for_iter(i) -> list(tree, tree_id, events[edge, is_SA, nchr, sys])

# dated trees and maps, with a hint on how to make them if missing
.TREES <- "tree/odon.boot.dated.tre"
.MAPS  <- "simmaps.rds"
if (!file.exists(.TREES))
  stop(paste0("tree/odon.boot.dated.tre not found.\n",
    "  It comes from 09_date.R (which needs support/ from 08b_support.sbatch).\n",
    "  Run:  sbatch 09_date.sbatch"), call. = FALSE)
if (!file.exists(.MAPS))
  stop(paste0("simmaps.rds not found. The stochastic maps come from:\n",
    "    Rscript 09b_prep_simmap.R\n",
    "    Rscript 09d_make_todo.R <n_total>\n",
    "    sbatch --array=1-$(wc -l < todo.txt)%200 09c_simmap.sbatch\n",
    "    Rscript 09e_collect.R"), call. = FALSE)

.trees <- ape::read.tree(.TREES)
.maps  <- readRDS(.MAPS)

# states "<nchr>_X0" / "<nchr>_neo" (09b_prep_simmap.R) before and after a transition ->
# list(is_SA, nchr before, sys before) for a fusion, NULL for fissions and reversals
classify_transition <- function(from, to) {
  f <- strsplit(as.character(from), "_", fixed = TRUE)[[1]]
  t <- strsplit(as.character(to),   "_", fixed = TRUE)[[1]]
  if (length(f) < 2 || length(t) < 2) return(NULL)
  n_from <- suppressWarnings(as.numeric(f[1]))
  n_to   <- suppressWarnings(as.numeric(t[1]))
  if (is.na(n_from) || is.na(n_to) || n_to >= n_from) return(NULL)  # not a fusion
  list(is_SA = as.integer(f[2] != t[2]), nchr = n_from,
       sys   = if (f[2] == "X0") "XO" else "XY")
}

# fusion events of one map; $maps[[e]] (same order as tr$edge rows) holds the states along edge e
extract_events <- function(m) {
  stopifnot(!is.null(m$maps))
  out <- list()
  for (e in seq_along(m$maps)) {
    seg <- names(m$maps[[e]])
    if (length(seg) < 2) next
    for (k in seq_len(length(seg) - 1L)) {
      cl <- classify_transition(seg[k], seg[k + 1L])
      if (is.null(cl)) next
      out[[length(out) + 1L]] <- data.table::data.table(
        edge = e, is_SA = as.integer(cl$is_SA), nchr = as.numeric(cl$nchr),
        sys = as.character(cl$sys))
    }
  }
  if (!length(out)) return(data.table::data.table(
    edge = integer(0), is_SA = integer(0), nchr = numeric(0), sys = character(0)))
  data.table::rbindlist(out)
}

n_iter_available <- function() length(.maps)

# map i with its tree; if the dated tree's edge matrix differs from the map's, use the map itself
events_for_iter <- function(i) {
  m   <- .maps[[i]]
  tid <- ((i - 1L) %% length(.trees)) + 1L
  tr  <- .trees[[tid]]
  if (!identical(tr$edge, m$edge)) {
    tr <- m; tid <- -i
  }
  list(tree = tr, tree_id = tid, events = extract_events(m))
}
