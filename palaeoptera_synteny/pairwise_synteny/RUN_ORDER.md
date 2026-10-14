# Figure 1 pipeline — run order

Everything runs from **`/path/to/project/palaeoptera_synteny/pairwise_synteny`**. Each script has a
`PATHS` block at the top with absolute paths; edit there, nowhere else.
The R scripts `setwd(WORKDIR)` themselves, so it doesn't matter where you
launch them from.

Nothing computed before steps 0–1 is valid afterwards: chromosome selection
changes `na`/`nb` and therefore every metric; the new tree changes `dist`,
the shared-range bound, and the residual strata.

| step | command | needs | writes |
|---|---|---|---|
| 0 | `python3 00_chromosomes.py` | `datasets` CLI, the existing download dir, `bed/` | `assemblies.tsv`, `chromosomes.tsv`, `seqrep/` |
| 1 | `sbatch 01_build_tree.sbatch` | OrthoFinder results, MAFFT, IQ-TREE | `tree/iqtree/tree.treefile` |
| 2 | `Rscript 02_synteny_pairs.R o2o`<br>`Rscript 02_synteny_pairs.R universal` | 0, 1, OrthoFinder results | `pairs_o2o.tsv` (primary), `pairs_universal.tsv` |
| 3 | `Rscript 03_stats.R` | `pairs_o2o.tsv` | console: every number for the Results |
| 3b | `Rscript 03b_compare_marker_sets.R` | both pair tables | `marker_set_comparison.tsv` (supplement) |
| 4 | `Rscript fig1bc.R` (in `figs/fig1/`, laptop) | `pairs_o2o.tsv` | `fig1b_retained.svg`, `fig1c_ari.svg` |
| 4b | `Rscript fig1a_tree.R` (in `figs/fig1/`, laptop) | `tree/iqtree/tree.treefile` | `palaeoptera.tree.svg` |

Step 0 reuses the metadata from your `download.sbatch` run
(`raw/ncbi_dataset/data/assembly_data_report.jsonl`) to map each genome to its
GCA accession, then pulls per-assembly sequence reports — a small metadata
download, not the genomes again. Check `assemblies.tsv` before continuing: any
genome it warns about is one where the bed basename didn't match an organism
name.

Step 1 is the only long job — submit it with `sbatch 01_build_tree.sbatch`
rather than running it on the login node. Check the allocation, core count and
walltime in the header first. Alignments are cached in `tree/aln/`, so a
resubmitted job re-uses them; `--resume` in the job script lets IQ-TREE pick up
from its checkpoint if the job hits walltime.

## Script output that belongs in the text
- `03_stats.R` §2b — what the excluded deep mayfly pairs would add, i.e. that the bound is conservative
- `03_stats.R` §4 — species-level non-overlap for displacement; overlap for ARI, and why both metrics are shown
- `03_stats.R` §5 — whether mayfly mixing is genome-wide or confined to a few chromosomes
- `03_stats.R` §7 — the three standing caveats
