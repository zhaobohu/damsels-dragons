# Panel B: neo-X / neo-Y divergence dating

Every cluster command runs from `/path/to/project/divergence_dating`. Steps 01–12 run in
order. `mcmctree.sh` is a helper called by 08, 10 and 12; never run it by hand.

## 0. Status (27 Sep 2026)

- **01–08 are done.** Rebuilding the locus selection from `gene_pairs_ds`
  (rank-filtered anchors, `min.ratio = 0.65`) reproduced the same 32 loci and
  byte-identical alignments (`ajun_concat.phy` md5 41571881e563ba8d5e3ae7d397f8e9ac).
  Full data, chain 1: t_n15 mean 2.97 Ma, 95% HPD 1.67–4.46 Ma.
- **The controls are being expanded from 99 to 999 per subset.** Only 09, 10
  and 11 rerun, then the plot. Section 4 has the commands.

## 1. Files

| file | role | run by | changed in this update |
|---|---|---|---|
| `01.cds.sbatch` | longest-isoform CDS for 10 species (AGAT) | `sbatch` | |
| `02.select_loci.R` | candidate loci: QC'd X–Y pairs ∩ rank-filtered anchors ∩ neo-X arm ∩ 1:1 in 9 outgroups | `Rscript` | |
| `anchor_filter.R` | sourced by 02: anchor-partner rank statistic, shared with 4a_new.R and fate.R | — | |
| `03.extract_loci.R` | one 11-tip FASTA per locus | `Rscript` | |
| `04.align_trees.sbatch` | PRANK codon alignment + IQ-TREE gene trees (seed 67) | `sbatch` | |
| `05.filter_loci.R` | sister tips, clean neo-Y CDS, ≥ 100 codons → `locus_table_final.tsv` | `Rscript` | |
| `06.build_phylip.R` | codon-partitioned concatenation + top 75/50/25% dS subsets | `Rscript` | |
| `07.modelsel.sh` | ModelFinder restricted to BASEML models; exits 1 unless GTR at all 3 positions | `bash` | |
| `08.inferences.sbatch` | six MCMCtree runs (2 full chains, 3 subsets, prior only), GTR+Γ5 | `sbatch` | |
| `09.random_controls.R` | **999** size-matched random subsets per fraction (seed 42) | `Rscript` | **yes** |
| `10.controls.sbatch` | 2997 control runs packed into a 500-task array, ≤ 45 running | `sbatch` | **yes** |
| `11.collect.sh` | ages table + permutation p per fraction | `bash` | **yes** |
| `12.sensitivity.sh` | optional HKY85+Γ5 comparison | `sbatch` | **yes** |
| `mcmctree.sh` | helper: per-dataset Hessian (usedata 3), then approximate likelihood (usedata 2) | — | |
| `mcmctree.ctl` | control-file template (`model = 7`) | — | |
| `ages.awk` | parses `t_n*` lines from `out.txt` | — | |
| `4b_forest.R` | the forest plot, run on the laptop | `Rscript` | **yes** |

Also needed: `ajun_calib.tre`, the pcm conda env, the genomes and GFFs under
`/path/to/project/data/`, the completed `gene_pairs_ds/`, and the
GENESPACE pangenes table.

## 2. What changed in this update

- **09.random_controls.R**
  - `n_rep` is now 999. The seed is still 42, so `rand75_001`–`099` are the
    same draws as before. Every other control is new, because draws are
    sequential.
  - Refuses to run while any old control files are present, so two control
    sets can't mix.
  - Stops if `n_rep > 999`, because the `rand??_???` IDs only have three digits.
  - Reports how many distinct locus sets each fraction drew.
- **10.controls.sbatch**
  - A fixed 500-task array (`--array=1-500%45`). Task *t* runs manifest rows
    *t*, *t*+500, *t*+1000, …, which is six controls for most tasks and two
    from each fraction.
  - The six run side by side on six CPUs, so a task takes about one control's
    time.
  - Each control writes its own log to `logs/controls/<id>.log`.
  - A task exits non-zero and names any control that failed.
  - **Safe to resubmit**: a control whose `out.txt` is complete and newer than
    its `.phy` is skipped. `FORCE=1` reruns everything.
- **11.collect.sh**
  - Adds an `n` column.
  - Warns if fewer controls were read than the manifest lists, because
    otherwise n and p would shrink without notice.
  - Prints rows in the fixed order 75 / 50 / 25.
  - No longer runs with `set -x`, which would print a trace over ~3000
    directories.
- **12.sensitivity.sh**: the two HKY chains now run in parallel, and the
  script has `#SBATCH` headers. Submit it with `sbatch`; don't run it on the
  login node.
- **4b_forest.R**
  - Control points drawn at alpha 0.08 with fixed vertical jitter, so 999
    points per lane show a density rather than a bar.
  - The header comment now says that point estimates are means.
  - Messages name `out.txt` rather than `mcmc.txt` when `STAT = "mean"`.

## 3. Full run order (reference; 01–08 are already done)

After each step: `echo "$(date -Is)  NN  <one-line result>" >> RUNLOG.md`.

```bash
sbatch 01.cds.sbatch                          # grep -c '^>' *.longest.cds.fna: 10 files, none empty
Rscript 02.select_loci.R                      # "158 -> 153 after rank filter"; -> locus_table.tsv
Rscript 03.extract_loci.R                     # -> per_locus/<pgID>/<pgID>.fasta
sbatch 04.align_trees.sbatch                  # .best.fas and .treefile counts must match
Rscript 05.filter_loci.R                      # wc -l < locus_table_final.tsv  -> 33 (32 loci)
Rscript 06.build_phylip.R                     # full 17741 | top75 10219 | top50 7010 | top25 2357 nt
bash 07.modelsel.sh                           # gate: GTR at all three positions
sbatch 08.inferences.sbatch                   # column -t modelsel/ages_all_runs.tsv
Rscript 09.random_controls.R                  # 2997 alignments + manifest
sbatch 10.controls.sbatch                     # 2997 control runs
bash 11.collect.sh                            # ages table + permutation p
sbatch 12.sensitivity.sh                      # optional
```

The neo-X/neo-Y node is `t_n15`. This was confirmed from the FigTree block of
`full_chain1/out.txt` and is hard-coded in 11, 12 and `4b_forest.R`.

## 4. Expanding the controls to 999 (from the current directory)

**Step 1: check timing and limits (one minute).**
```bash
cd /path/to/project/divergence_dating
# how long one control took last time (elapsed per array task, longest last)
sacct -u "$USER" --name=controls -S 2026-09-01 -X --format=JobID%20,Elapsed,State | sort -k2 | tail -3
scontrol show config | grep -i MaxArraySize   # must be > 500
```
Expected wall time is about 12 × the longest elapsed above, plus queue wait
(section 5).

**Step 2: move the 99-control set aside.**
```bash
mkdir -p controls_n99
mv rand??_??? ajun_concat_rand*.phy controls_manifest.tsv controls_n99/
cp modelsel/ages_controls.tsv controls_n99/
echo "$(date -Is)  ctl  99-control set moved to controls_n99/" >> RUNLOG.md
```

**Step 3: install the updated scripts** (`4b_forest.R` goes to the laptop).
```bash
cp /path/to/update/{09.random_controls.R,10.controls.sbatch,11.collect.sh,12.sensitivity.sh} .
md5sum *.R *.sh *.sbatch *.awk *.ctl > checksums.md5
echo "$(date -Is)  ctl  scripts updated for 999 controls; checksums.md5 regenerated" >> RUNLOG.md
```

**Step 4: generate the alignments** (a few minutes, on the login node).
```bash
Rscript 09.random_controls.R                  # "999 controls written, N distinct locus sets" x 3
wc -l < controls_manifest.tsv                 # 2998 = 2997 + header
```

**Step 5: submit.** There is one array job of 500 tasks, with at most 45
running at once.
```bash
sbatch 10.controls.sbatch
echo "$(date -Is)  10  submitted job <JOBID>" >> RUNLOG.md
```
Monitor:
```bash
squeue -u "$USER" --name=controls -h -t R  | wc -l                  # running tasks (<= 45)
squeue -u "$USER" --name=controls -h -t PD -o '%r' | sort | uniq -c  # why tasks are pending
grep -l '^t_n15' rand??_???/out.txt 2>/dev/null | wc -l             # finished controls, of 2997
```

**Step 6: when the array has finished, check for failures and resubmit if needed.**
Don't resubmit while the first array is still running. A control that is
still in progress has no finished `out.txt` yet, so it would be started a
second time in the same directory.
```bash
sacct -j <JOBID> -X --format=JobID%20,State,Elapsed | grep -v COMPLETED
cat logs/controls.<JOBID>_*.err | grep FAILED                       # which controls, if any
sbatch 10.controls.sbatch                                           # reruns only what is missing
```
Finished controls are skipped, so a resubmitted task with nothing left to do
exits in seconds. To rerun only specific tasks: `sbatch --array=41,203 10.controls.sbatch`.

**Step 7: collect, and check against the old set.**
```bash
bash 11.collect.sh                           # n must be 999 in all three rows; no WARNING
# rand75_001-099 are the same alignments as before. If mcmctree.ctl is unchanged
# since the 99-control run, their results should be identical (no output here):
for d in controls_n99/rand75_0??; do
  b=$(basename "$d")
  cmp -s <(grep '^t_n15' "$d/out.txt") <(grep '^t_n15' "$b/out.txt") || echo "differs: $b"
done | head
echo "$(date -Is)  11  999 controls collected: <paste the three p values>" >> RUNLOG.md
```

**Step 8: plot (laptop).** Sync only what `4b_forest.R` reads, which is each
run's `out.txt` (not the Hessian's) and the manifest:
```bash
rsync -av --prune-empty-dirs \
  --exclude='hessian/' --exclude='per_locus/' --exclude='logs/' --exclude='controls_n99/' \
  --include='*/' --include='out.txt' --include='controls_manifest.tsv' --exclude='*' \
  <user>@<cluster>:/path/to/project/divergence_dating/ <laptop>/divergence_dating/
Rscript 4b_forest.R                           # prints the [controls] table, writes the svg
```
The `[controls]` table printed by `4b_forest.R` should match 11's output.

## 5. Resources, and how to change them

| setting | value | why |
|---|---|---|
| `--array` | `1-500%45` | 500 tasks, never more than 45 running |
| `STRIDE` (in the script) | 500 | must equal the array size; rows per task = ceil(2997 / STRIDE) = 6 |
| `--cpus-per-task` | 6 | one CPU per control, so a task's six controls run at once |
| `--time` | 02:00:00 | covers one control per CPU; the same limit a single control had before |
| `--mem-per-cpu` | 4G | unchanged; MCMCtree needs far less |

- **Wall time.** 500 tasks at 45 concurrent is ceil(500/45) = 12 rounds, and
  each round takes about one control's time, so the total is about 12 × one
  control plus queue wait. For comparison, the 297-control run needed 15
  rounds at 20 concurrent, one control each. The 999 set should therefore
  finish in about the same wall time as the old 297, if the queue gives you
  the cores.
- **Cores.** 45 tasks × 6 CPUs = 270 cores at peak. If step 5 shows tasks
  pending with a reason like `QOSMaxCpuPerUserLimit` or `AssocGrpCpuLimit`,
  cancel, set `--cpus-per-task=3` and `--time=04:00:00` (two rounds of
  controls per task), and resubmit. Finished controls are skipped.
- **Task limit.** The array is exactly 500 tasks. If your 500 limit counts
  every job you have queued, other pending jobs count against it too.
- **To use a different array size**, change `--array=1-N%45` and `STRIDE=N`
  together. Then set `--cpus-per-task` to ceil(2997 / N).

## 6. Manuscript text to update after step 7

- **Methods, random controls paragraph:**
  - “we generated for each dS subset 99 random control subsets” → 999
  - “where n = 99 controls” → n = 999
  - Optionally add “(seed 42)”.
- **Fig. 4 legend:** “Grey points beneath each subset are the 99 size-matched
  random control analyses” → 999.
- **Results, from 11's output:**
  - control means (now 3.63 / 3.62 / 4.41 Ma)
  - separations (now 2.27 / 2.85 SD; 1.66 for the top 25%) = dS_eff / ctl_sd
  - p-values (now 0.020 / 0.010 / 0.080)
  - the control SDs in the model-selection paragraph (now 1.18 vs 0.67 Ma)
- The smallest attainable p is now 0.001. If a subset has k = 0, report
  p = 0.001 and say no control reached the observed age.

## 7. Settings in force (for the methods)

- PRANK `-codon -F`, seed 67; IQ-TREE 3 `-m MFP -B 1000`, seed 67
- locus filters: neo-X and neo-Y sister in the unrooted ML gene tree; neo-Y CDS
  length divisible by 3 with no internal stop; ≥ 100 codons
- MCMCtree:
  - three codon-position partitions; GTR+Γ5 (`model = 7`, `alpha = 0.5`, `ncatG = 5`)
  - independent-rates clock (`clock = 2`)
  - priors: `rgene_gamma = 2 20 1`, `sigma2_gamma = 1 10 1`, `BDparas = 1 1 0.1`
  - burn-in 200,000, then 20,000 samples every 100 iterations
  - per-dataset Hessian
  - seeds: 67 and 69 for the two full chains; 67 for everything else
- controls:
  - 999 per fraction
  - locus count matched; total length within ± 15%
  - seed 42
  - p = (1 + k) / (n + 1)

## 8. Provenance

`checksums.md5` (regenerated in step 3) records which script versions ran.
`mcmctree.sh` writes `provenance.md5` into every run directory (alignment md5,
Hessian md5, the ctl lines that matter). `RUNLOG.md` records what was run and
when.
