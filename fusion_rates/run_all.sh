#!/usr/bin/env bash
# The full pipeline as a checklist, run from inside the run directory. Several steps need their
# output inspected before continuing. For a new tree build use the driver instead:
#   bash setup_new_run.sh <old_run_dir> <new_run_dir>
#   cd <new_run_dir> && tmux new -s tree && bash run_tree.sh
set -euo pipefail

# names
Rscript 00_resolve_names.R          # GBIF -> name_map.csv
#   INSPECT the FUZZY block and the genus reassignments before continuing.

# sequences (01 rewrites prep.csv and unfiltered.COI.fna, so always rerun 02 after it; 02 is cached)
Rscript 01_prep.R                   # -> taxa.tsv, unfiltered.COI.fna, prep.csv
Rscript 02a_seed_cache.R            # ONCE: import already-downloaded sequences
Rscript 02_fetch.R                  # NCBI: COI gaps + non-COI markers
#   If every chunk reports 0 records, the search is failing, not the data.
Rscript 03_screen.R                 # unambiguous length + stop-free frame
sbatch  04_blast.sbatch             # per-marker reference + BLAST coverage
#   wait for hits_<marker>.tsv; check refs/reference_report.tsv

Rscript 05_select.R                 # coverage filter -> <marker>.fna
#   INSPECT scov.pdf and coverage_report.tsv before accepting SC_MIN.

Rscript 05b_subset.R                # restrict to MARKERS_TREE (config.R)

# alignment and constraint
bash    06_align_trim.sh
#   CHECK align/occupancy.tsv -- drop any marker below ~0.6
Rscript 07_constraint.R
Rscript 07b_check.R                 # validate before spending hours on 08
#   MUST exit 0: constraint/alignment labels, charset ranges, neo-sex tips in the focal clade.

# tree, support and dating
sbatch  08_tree.sbatch              # 10-seed array
#   wait, then pick the best topology:
#     grep -H "BEST SCORE" run.*/run.*.log | sort -t: -k3 -g
sbatch  08b_support.sbatch          # uses the best-scoring run.N unless BEST_TREE is exported
#   CHECK: grep -c ';' support/support.ufboot  -- must be 2000, not 0

sbatch  09_date.sbatch              # chronos on every bootstrap tree; NOT on a login node
#   Checkpoints to tree/boot_dated_ckpt.rds every chunk -- a timeout resumes.

# stochastic maps
Rscript 09b_prep_simmap.R           # states, rate matrix, one Q fit -> prep.rds
#   NOTE which fit it used. If it fell back to a flat prior for the fit, report that.
Rscript 09d_make_todo.R "$N_TOTAL"   # N_TOTAL is printed by 09b
sbatch --array=1-$(wc -l < todo.txt)%200 09c_simmap.sbatch
#   ONE map per task. Rerun 09d_make_todo.R after the array ends and resubmit
#   for any that failed.
Rscript 09e_collect.R               # maps/ -> simmaps.rds
#   CHECK the SA count per map is not 0 (09e warns if it is).

# time-slice test, Fig 2C table and the remaining numbers
sbatch  15_timeslice.sbatch         # time-slice test: all maps x 6 T factors x 9,999 relabellings
#   -> timeslice_clades.tsv, timeslice_permap.tsv, timeslice_summary.tsv
#   (read by figs/fig2/fig2d.R and supp_cut_age.R). Directly: Rscript 15_timeslice.R [n_maps] [n_null] [cores]
sbatch  13_figures.sbatch           # after 15_timeslice: --dependency=afterok:<its job id>
#   13_figures.R on all maps -> fig_levels_permap.tsv, read by figs/fig2/fig2c.R;
#   then get_numbers.R clades -> remaining numbers for the text, in the job log.
