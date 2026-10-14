#!/usr/bin/env Rscript
# One stochastic map per SLURM array task (map index = line SLURM_ARRAY_TASK_ID of todo.txt)

source("config.R")
source("simmap_uniform.R")

# map index for this task; exit quietly past the end of todo.txt or if the map exists
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
if (is.na(task_id)) stop("No SLURM_ARRAY_TASK_ID -- this is an array task script.")

if (!file.exists("todo.txt")) stop("todo.txt not found -- run 09d_make_todo.R first.")
todo <- as.integer(readLines("todo.txt"))
if (task_id > length(todo)) {
  message(sprintf("task %d: past the end of todo.txt (%d entries) -- nothing to do.",
                  task_id, length(todo)))
  quit(save = "no", status = 0)
}
g <- todo[task_id]

dir.create("maps", showWarnings = FALSE)
out <- sprintf("maps/simmap_part_%04d_%04d.rds", g, g)
if (file.exists(out)) {
  message(sprintf("task %d: map %d already on disk -- skipping.", task_id, g))
  quit(save = "no", status = 0)
}

# tree for this map (NSIM_PER_TREE maps per tree)
prep     <- readRDS("prep.rds")
tree_idx <- ((g - 1L) %/% NSIM_PER_TREE) + 1L
if (tree_idx > length(prep$trees))
  stop(sprintf("map %d wants tree %d but prep.rds holds %d",
               g, tree_idx, length(prep$trees)))
message(sprintf("task %d: map %d (tree %d), seed %d", task_id, g, tree_idx, SEED0 + g))

# seeded by map index, so map g is reproducible on its own
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED0 + g)
t0 <- Sys.time()

# one map, endpoint-conditioned by uniformization (simmap_uniform.R)
tr <- prep$trees[[tree_idx]]
xx <- prep$x
s <- simmap_uniform(tr, xx, prep$Qfit, prep$pi_vec, nsim = 1)
if (is.null(s[["maps"]])) stop(sprintf("map %d: result has no $maps", g))

# write to a per-task temp name, then rename into place
tmpf <- sprintf("%s.%s_%s.part", out,
                Sys.getenv("SLURM_ARRAY_JOB_ID", "na"),
                Sys.getenv("SLURM_ARRAY_TASK_ID", "na"))
saveRDS(list(s), tmpf, compress = "xz")
file.rename(tmpf, out)
message(sprintf("task %d: wrote %s in %.2f h", task_id, out,
                as.numeric(difftime(Sys.time(), t0, units = "hours"))))
