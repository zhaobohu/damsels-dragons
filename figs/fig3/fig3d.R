library(data.table)

# neo-X / neo-Y divergence date (MCMCtree posterior mean, 95% HPD) for each dS subset,
# against its size-matched random controls

RUNS  <- "../../divergence_dating"
AGES  <- file.path(RUNS, "modelsel", "ages_controls.tsv")   # from 11.collect.sh
MAN   <- file.path(RUNS, "controls_manifest.tsv")            # from 09.random_controls.R
NODE  <- "t_n15"                                  # (neo-X, neo-Y)
SCALE <- 100                                      # time units (100 Ma) -> Ma

subsets <- data.frame(
  key   = c("full",           "rand75",            "rand50",            "rand25"),
  dir   = c("full_chain1",    "ds_top75",          "ds_top50",          "ds_top25"),
  label = c("Full (32 loci)", "Top 75% (24 loci)", "Top 50% (16 loci)", "Top 25% (8 loci)"),
  stringsAsFactors = FALSE)
ctl_keys <- c("rand75", "rand50", "rand25")       # subsets with random controls

OUT       <- "neoXY_date_forest.svg"
col_ctl   <- "grey45"
col_ctlmn <- "grey20"
ctl_off   <- 0.22        # control lane sits this far below its row
ctl_alpha <- 0.1

# node ages for every run, in Ma
ages <- fread(AGES)[node == NODE]
ages[, `:=`(point = mean * SCALE, hpd.lower = hpd_lo * SCALE, hpd.higher = hpd_hi * SCALE)]

# observed runs
obs <- ages[match(subsets$dir, run)]
if (anyNA(obs$run))
  stop("not in ", AGES, ": ", paste(subsets$dir[is.na(obs$run)], collapse = ", "),
       "\n  Rerun divergence_dating/11.collect.sh after these runs finish.")
dat <- as.matrix(obs[, .(point, hpd.lower, hpd.higher)])
rownames(dat) <- subsets$key

cat("\n[observed] posterior mean and 95% HPD (Ma)\n")
print(round(dat, 3))

# control ages, matched to the manifest
man <- fread(MAN)
man[, subset := sub("_[0-9]+$", "", id)]
man[ages, age := i.point, on = .(id = run)]

n_missing <- man[is.na(age), .N]
if (n_missing == nrow(man))
  stop("no control runs (", man$id[1], " ...) found in ", AGES)
if (n_missing)
  warning(n_missing, " of ", nrow(man), " control run(s) listed in ", MAN,
          " are missing from ", AGES, " -- dropped from the lanes and the statistics below")
message(sprintf("[controls] read %d of %d control runs", nrow(man) - n_missing, nrow(man)))

# control statistics; p = (1 + k) / (n + 1), k = controls with age >= observed (ties included)
cat("\n[controls] posterior means (Ma)\n")
ctl_stats <- list()
for (k in ctl_keys) {
  a     <- man[subset == k & !is.na(age), age]
  o     <- dat[k, "point"]
  n_ge  <- sum(a >= o)
  ctl_stats[[k]] <- data.table(subset = k, n = length(a),
                               ctl_mean = mean(a), ctl_sd = sd(a), ctl_median = median(a),
                               observed = o,
                               sd_above = (o - mean(a)) / sd(a),
                               n_ge_obs = n_ge, n_tie = sum(a == o),
                               p_perm = (1 + n_ge) / (length(a) + 1))
}
ctl_stats <- rbindlist(ctl_stats)
print(ctl_stats[, lapply(.SD, function(z) if (is.numeric(z)) round(z, 3) else z)])

# plot layout
n       <- nrow(dat)
y.pos   <- n:1                                   # Full at the top
val.lab <- sprintf("%.2f [%.2f, %.2f]",
                   dat[, "point"], dat[, "hpd.lower"], dat[, "hpd.higher"])
x.max   <- max(dat[, "hpd.higher"], man$age, na.rm = TRUE) * 1.12

y.ctl <- setNames(y.pos[match(ctl_keys, subsets$key)], ctl_keys)

svg(OUT, width = 28, height = 8, bg = "transparent")
par(mar = c(4, 11, 2, 1))
plot(NA, xlim = c(0, x.max), ylim = c(0.5, n + 0.5),
     xlab = "neo-X / neo-Y divergence (Ma)", ylab = "", yaxt = "n", bty = "n",
     mgp = c(3, 3.5, 0), tcl = -2)

# controls first, in their own lane below each row: point estimates from separate runs
# (a null for dataset size), not an interval on the observed date
for (k in names(y.ctl)) {
  a <- man[subset == k & !is.na(age), age]
  if (!length(a)) next
  points(a, rep(y.ctl[[k]] - ctl_off, length(a)),
         pch = 16, cex = 1.6, col = adjustcolor(col_ctl, ctl_alpha))
  segments(mean(a), y.ctl[[k]] - ctl_off - 0.10,
           mean(a), y.ctl[[k]] - ctl_off + 0.10, col = col_ctlmn, lwd = 4)
}

# observed posterior means with 95% HPD bars
segments(dat[, "hpd.lower"], y.pos, dat[, "hpd.higher"], y.pos, lwd = 9)
segments(dat[, "hpd.lower"],  y.pos - 0.12, dat[, "hpd.lower"],  y.pos + 0.12, lwd = 9)
segments(dat[, "hpd.higher"], y.pos - 0.12, dat[, "hpd.higher"], y.pos + 0.12, lwd = 9)
points(dat[, "point"], y.pos, pch = 19, cex = 5)

axis(2, at = y.pos, labels = subsets$label, las = 1, tick = FALSE)
text(dat[, "hpd.higher"], y.pos, val.lab, pos = 4, cex = 0.75)

dev.off()
message("[wrote] ", OUT)
