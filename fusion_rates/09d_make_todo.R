#!/usr/bin/env Rscript
# Write todo.txt, the map indices not yet in maps/. Run only after the previous array has fully ended.
#   Rscript 09d_make_todo.R <n_total>

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("usage: Rscript 09d_make_todo.R <n_total>")
n_total <- as.integer(args[1])

f <- list.files("maps", pattern = "^simmap_part_[0-9]+_[0-9]+\\.rds$")
if (!length(f)) message("maps/ is empty -- todo will be every index.")

# mark every index covered by a part file (simmap_part_<first>_<last>.rds)
have <- logical(n_total)
m <- regmatches(f, regexec("^simmap_part_([0-9]+)_([0-9]+)\\.rds$", f))
for (v in m) {
  a <- as.integer(v[2]); b <- as.integer(v[3])
  if (a < 1L || b > n_total || a > b)
    stop(sprintf("part file with a nonsensical range: %d..%d", a, b))
  have[a:b] <- TRUE
}

todo <- which(!have)
writeLines(as.character(todo), "todo.txt")
message(sprintf("%d part files -> %d/%d maps present, %d missing.",
                length(f), sum(have), n_total, length(todo)))

# report the missing indices as runs ("3", "5-9", ...)
if (length(todo)) {
  brk <- c(0, which(diff(todo) != 1L), length(todo))
  rng <- character(0)
  for (i in seq_len(length(brk) - 1L)) {
    s <- todo[brk[i] + 1L]; e <- todo[brk[i + 1L]]
    rng <- c(rng, if (s == e) as.character(s) else sprintf("%d-%d", s, e))
  }
  message("missing: ", paste(utils::head(rng, 20), collapse = ", "))
  message(sprintf("\nsubmit with:  sbatch --array=1-%d%%200 09c_simmap.sbatch", length(todo)))
} else message("nothing missing -- run 09e_collect.R")
