#!/usr/bin/env Rscript
# Build prep.rds (trees, tip states, fitted Q, root prior) for the per-map array.
# States are "<haploid nchr>_<X0|neo>", as parsed by classify_transition() in events_adapter.R.

suppressPackageStartupMessages({
  library(ape); library(data.table); library(phytools)
})
source("config.R")

TREES <- "tree/odon.boot.dated.tre"
if (!file.exists(TREES)) stop(TREES, " not found -- run 09_date.R first", call. = FALSE)

tx <- fread(TAXA_TSV)
trees <- read.tree(TREES)
if (inherits(trees, "phylo")) trees <- list(trees)

# tip states for karyotyped tips in the tree
tx[, sexcls := fifelse(sex_system == "X", "X0", "neo")]
tx[, tip := gsub(" ", "_", species)]
tx <- tx[!is.na(nchr) & tip %chin% trees[[1]]$tip.label]

tx[, state := sprintf("%d_%s", as.integer(nchr), sexcls)]

# contiguous state space, unobserved intermediate counts included (keeps the rate graph connected)
rng <- seq(min(tx$nchr), max(tx$nchr))
lv  <- as.vector(t(outer(rng, c("X0", "neo"), paste, sep = "_")))
n_of <- as.integer(sub("_.*", "", lv)); s_of <- sub(".*_", "", lv)

x <- matrix(0, nrow(tx), length(lv), dimnames = list(tx$tip, lv))
x[cbind(seq_len(nrow(tx)), match(tx$state, lv))] <- 1

# rate classes: 1 = AA fusion, 2 = SA fusion (X0 -> neo), 3 = fission (one step each),
# 4 = SA-fusion reversal (neo -> X0 at unchanged count)
Q <- matrix(0L, length(lv), length(lv), dimnames = list(lv, lv))
for (i in seq_along(lv)) for (j in seq_along(lv)) {
  if (i == j) next
  dn <- n_of[i] - n_of[j]
  if (dn == 1L) {
    if (s_of[i] == "X0" && s_of[j] == "neo") Q[i, j] <- 2L
    else if (s_of[i] == s_of[j])             Q[i, j] <- 1L
  } else if (dn == -1L && s_of[i] == s_of[j]) Q[i, j] <- 3L
  else if (dn == 0L && s_of[i] == "neo" && s_of[j] == "X0")
    Q[i, j] <- 4L
}
stopifnot(any(Q == 1L), any(Q == 2L))

cat(sprintf("%d tips, %d states over n = %d-%d (%d observed)\n",
            nrow(tx), length(lv), min(rng), max(rng), uniqueN(tx$state)))
cat(sprintf("model: %d AA, %d SA, %d fission, %d reversal\n",
            sum(Q == 1L), sum(Q == 2L), sum(Q == 3L), sum(Q == 4L)))

# ML tree (or the first bootstrap tree) restricted to karyotyped tips
ml <- if (file.exists("tree/odon.dated.tre")) read.tree("tree/odon.dated.tre") else trees[[1]]
ml <- keep.tip(ml, intersect(ml$tip.label, rownames(x)))
cat("fitting Q once on the ML tree (reused for every map)...\n")

# root prior: X0 at every count, weighted by that count's frequency among X0 tips
x0_counts <- table(factor(tx[sexcls == "X0"]$nchr, levels = rng))
root_pi <- rep(0, length(lv)); names(root_pi) <- lv
root_pi[sprintf("%d_X0", rng)] <- as.numeric(x0_counts) / sum(x0_counts)
nz <- root_pi[root_pi > 0]
cat(sprintf("root prior (empirical_X0): %s\n",
            paste(sprintf("%s=%.2f", names(nz), nz), collapse = " ")))

# fit Q under the root prior; fall back to a flat prior if that fails
fit <- NULL; used <- NA_character_
for (att in list(list(p = root_pi, lab = "fusions + fissions, pi=empirical_X0"),
                 list(p = "equal", lab = "fusions + fissions, pi=equal"))) {
  fit <- tryCatch(fitMk(ml, x[ml$tip.label, , drop = FALSE], model = Q, pi = att$p),
                  error = function(e) e)
  if (inherits(fit, "error") || any(!is.finite(fit$rates))) fit <- NULL
  if (!is.null(fit)) { used <- att$lab; fit_pi <- att$p; break }
  cat(sprintf("  %s: failed, trying next\n", att$lab))
}
if (is.null(fit)) stop("no model configuration converged -- inspect the state space",
                       call. = FALSE)
cat(sprintf("  model used: %s\n", used))
if (grepl("pi=equal", used))
  cat("  NOTE: Q was FITTED under a flat prior (fallback). The maps still use\n",
      "  the empirical X0 root prior; only the rate estimates saw the flat prior.\n", sep = "")
cat(sprintf("  rates: qAA = %.4g, qSA = %.4g%s\n", fit$rates[1], fit$rates[2],
            if (length(fit$rates) > 2) sprintf(", fission = %.4g", fit$rates[3]) else ""))
cat(sprintf("  qSA/qAA = %.3f; random pairing at n = %d would give %.3f\n",
            fit$rates[2] / fit$rates[1], round(median(tx$nchr)),
            (2 / median(tx$nchr)) / (1 - 2 / median(tx$nchr))))

# every tree restricted to karyotyped tips
kept <- vector("list", length(trees))
for (k in seq_along(trees)) kept[[k]] <- keep.tip(trees[[k]], intersect(trees[[k]]$tip.label, rownames(x)))
trees <- kept

# fitted rate matrix: model classes index fit$rates, diagonal = minus the row sum
Qnum <- matrix(0, length(lv), length(lv), dimnames = list(lv, lv))
nz <- which(Q > 0L, arr.ind = TRUE)
Qnum[nz] <- fit$rates[Q[nz]]
diag(Qnum) <- -rowSums(Qnum)
stopifnot(all(abs(rowSums(Qnum)) < 1e-8), all(Qnum[upper.tri(Qnum) | lower.tri(Qnum)] >= 0))

# maps always use the empirical X0 prior (pi_vec), whatever prior the fit needed (fit_pi)
pi_vec <- root_pi
if (!identical(fit_pi, root_pi))
  cat(sprintf("note: Q was fitted under pi=%s but maps use the empirical X0 prior -- intended\n",
              if (is.character(fit_pi)) fit_pi else "custom vector"))
saveRDS(list(trees = trees, x = x, tm = Q, pi = root_pi, fit_pi = fit_pi, pi_vec = pi_vec, Qfit = Qnum,
             levels = lv, model_used = used),
        "prep.rds")

n_total <- length(trees) * NSIM_PER_TREE
cat(sprintf("\nprep.rds written: %d trees x %d maps = %d total\n",
            length(trees), NSIM_PER_TREE, n_total))
cat(sprintf("next:\n  Rscript 09d_make_todo.R %d\n", n_total))
cat("  sbatch --array=1-$(wc -l < todo.txt)%200 09c_simmap.sbatch\n")
