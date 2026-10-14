# TOGA2 run: Anax junius (reference) -> Aeshna juncea (query)

Goal: classify every reference gene's projection in the *A. juncea* male assembly as intact /
partially intact / lost / missing, with inactivating mutations, so the fate.R gene set can be
filtered to validated genes. TOGA2 v2.0.10 + make_lastz_chains v3.1.7 (both hillerlab, GitHub).

Reference = the annotated genome whose transcripts are projected (Anax). Query = the genome being
tested (A. juncea). This matches fate.R, which is keyed by focal Anax gene.

## Files

| file | where it runs | what it does |
|---|---|---|
| `config.sh` | sourced by all | paths, allocation, CPU/memory budget. **Edit first.** |
| `00_setup_env.sh` | login node | clone+build TOGA2 and make_lastz_chains (+ container) into `/path/to/project/programs`; rerun-safe |
| `01_prep_genomes.sh` | login node | strip `.1` from sequence names, 2bit, chrom.sizes, softmask check, reference GFF -> renamed GTF |
| `02_make_chains.slurm` | `sbatch` | whole make_lastz_chains run in one 64-CPU x 3500 MB job |
| `02b_check_chains.sh` | login node / `srun` | job state, task failures, chain file integrity, aligned bases per chromosome and of the reference CDS |
| `03_prepare_ref.sh` | login node | `toga2.py prepare-input`: BED12, isoforms, U12 introns, CESAR2 profiles |
| `04_toga2_run.slurm` | `sbatch` | whole `toga2.py run` in one 64-CPU x 3500 MB job |
| `05_join.R` | login node | TOGA2 status per Anax gene, split by where each copy lands (neo-Y / X / elsewhere) + juncea drop list |

## Run order

```bash
cd /path/to/project/gene_fate
vi config.sh                        # 0. paths, ALLOC, NCPU/NMEM
# one time only -- the conda env (nothing below installs or updates it):
module load mamba/1.5.8 && mamba create -y -p /path/to/project/conda/envs/toga -c conda-forge -c bioconda 'python>=3.12,<3.13' 'nextflow>=25.10.0' 'rust>=1.90.0' gffread git wget make ucsc-fatotwobit ucsc-twobitinfo
bash 00_setup_env.sh                # 1. first run ~30-60 min; later runs skip what is already built
bash 01_prep_genomes.sh             # 2. minutes; confirm both genomes report a sensible softmasked %
sbatch 02_make_chains.slurm         # 3. hours to a couple of days; watch: tail -f mlc.*.out
bash 02b_check_chains.sh            #    did it work? COMPLETED 0:0, every step has COMPLETED rows,
                                    #    both name-mismatch counts 0, reference CDS mostly aligned
bash 03_prepare_ref.sh              # 4. minutes to ~1 h
sbatch 04_toga2_run.slurm           # 5. hours to days; watch: tail -f toga2.*.out
Rscript 05_join.R /path/to/project/gene_fate /path/to/project/data/genes/Aeshna.juncea.longest.gff   # 6.
```

Then, one filter for fate.R: `juncea_broken_genes.txt` is a DROP list for the juncea rows of the
pangenes file and the MCScanX pairs. "Overlapping" means a juncea CDS overlaps a projected exon on the
same strand (not just the projection's span, which also covers genes inside its introns). A juncea
gene is on it only if nothing overlapping it is better
than L and the L comes from an orthologous (not paralog-labelled) projection of an Anax gene that
projects intact somewhere else -- so a bad Anax model or a paralog's alignment can't delete a real
gene. One exception (gametolog rule): for an Anax gene on the neo-X homolog (OZ349811) whose copy on
one sex chromosome is intact, a paralog-labelled L copy on the other sex chromosome counts as evidence,
because TOGA2 gives the ortholog label to the intact chain and the other copy is the gametolog.
`juncea_gene_status.tsv` records the evidence in `l_src` (orthologous) and `gam_src` (gametolog), and
`kept_reason` for L-best genes that were not dropped. Don't use `juncea_valid_genes.txt` as a keep-list: it would also drop the ~30% of juncea genes
with no overlapping projection (genes with no Anax counterpart, or in sequence the chains don't cover).
`anax_unverified_genes.txt` (Anax genes with no intact copy anywhere) is information only: their
juncea copies are never dropped. Don't exclude their orthogroups -- that would also remove genes truly
lost from juncea, i.e. empty the "fully deleted" class, and it hits the anc-X homolog harder (28%
unverified vs 16.5% on autosomes).

`toga2_out/join/anax_gene_status.tsv` is the direct check on the "degenerated" class. One row per
Anax gene; `Y_best` / `Y_n` / `Y_all` are the statuses of its copies on the neo-Y only (likewise `X_*`,
`other_*`), `genome_wide` is TOGA2's gene-level status (best copy anywhere). `Y_call`:

| Y_call | TOGA2 status of the best neo-Y copy | reading |
|---|---|---|
| `intact` | FI / I | an intact copy on the neo-Y -- but it may be a paralog: see `anax_projections.tsv` |
| `partially_intact` | PI | no inactivating mutation, but part of the gene not aligned: ambiguous |
| `inactivated` | UL / L | copy present with inactivating mutations: degenerated (pseudogene) |
| `missing` | M / PM | most of the gene not aligned on the Y copy (deletion or assembly gap): ambiguous |
| `no_Y_copy` | none | no chain with score >= 15000 places the gene on the neo-Y: deleted, or too diverged to align |

`X_call` is the same call for the gene's copy on the X -- the matched baseline: for genes still on the
neo-X it shows how often TOGA2 calls a retained copy lost. Copies TOGA2 labels paralogous get their
L/M calls reported as `PG`; 05 converts those back to L or M from `meta/transcript_meta.tsv`, with the
thresholds TOGA2 itself uses (intact fraction < 0.35, non-deleted fraction < 0.2, or a mutation in the
middle 80% -> L; otherwise M).

`toga2_out/join/anax_projections.tsv` lists every aligned projection, one row per (projection,
juncea sequence): Anax gene, region (Y / neoX / ancX / auto / other), status, TOGA2's paralog and retro
labels, projected vs reference exon count, the juncea gene(s) whose CDS it overlaps, coordinates.

fate.R reads it from `TOGA_DIR` (with `anax_gene_status.tsv` and `juncea_broken_genes.txt`) for two things:
- **TOGA2 presence rescue** (`USE_TOGA_RESCUE`): a home-arm or sister-arm count of 0 is replaced by the
  number of the og's focal Anax genes that TOGA2 projects intact (FI/I), orthology-labelled, and not as
  a processed pseudogene (not `#retro`, not an intronless copy of a multi-exon gene) on that arm. On
  the autosomal baseline it applies to ogs with no juncea copy at all. The script prints the category
  table without -> with the rescue either way.
- **Breakdown of the "degenerated" genes** that have an intact neo-Y call (annotation-only
  classification): same og (should not happen), different og (og split), not PASS/array, not in
  pangenes, no gene model (BRAKER miss) -- all rescued -- versus retrocopy and paralog only, where the
  degenerated call stands. Per-projection detail goes to `degenerated_Y_intact.tsv`.

## Key outputs (`toga2_out/`)

- `loss_summary.tsv` -- `level entry status`; level = PROJECTION / TRANSCRIPT / GENE.
  Status codes: FI fully intact, I intact, PI partially intact, UL uncertain loss, L lost,
  M missing, PG paralogous projection, PP processed pseudogene, N no data.
- `inactivating_mutations.tsv` -- stop codons, frameshifts, splice-site and deletion events per projection.
- `orthology_classification.tsv` -- `t_gene t_transcript q_gene q_transcript orthology_class` (one2one etc.).
- `query_annotation.bed` -- projections in the final annotation, on juncea coordinates;
  name = `<transcript>#<gene>#<chain id>` (prepare-input names reference transcripts `<transcript>#<gene>`).
- `meta/loss_summary_extended.tsv`, `meta/query_annotation.with_discarded_exons.bed` -- every aligned
  projection with its own status, before filtering. 05_join.R reads these, not the filtered files.
- `summary.txt`, `meta/` (rejection logs etc.).

## Resume

- make_lastz_chains: rerun `sbatch 02_make_chains.slurm` (`-resume` is set) or use its
  `--from chain_antirepeat|fill_chains|clean_chains` checkpoints (README in the repo).
- TOGA2: same `-o toga2_out` plus `--resume_from <step>`; steps are listed in `04_toga2_run.slurm`.

## Things I guessed -- verify before running

- **Per-task memory in 02**: LASTZ/axtChain/RepeatFiller are given 3 GB each (pipeline asks 24-50 GB).
  A task killed for memory is retried with 6 GB, 9 GB, ... (not needed in the first run: no retries).
- **TOGA2 alignment split in 04**: 1 slot each for the 5/7/10/15 GB CESAR bins, the rest (54 at
  64 x 3500M) for the 3 GB bin. If `toga2.*.out` shows the 3 GB bin done long before the others,
  shift slots toward the bigger bins in the `nfcfg` lines.
- **LASTZ parameters**: left at the pipeline defaults (mammal-tuned). Anax and Aeshna are both
  Aeshnidae, so this should be fine; if chains look sparse (small `07_final` file, many M in
  loss_summary), loosen `lastz_y/l/k` in a copy of `params.json`.
- **Anax gene ids for the fate.R join**: TOGA2's gene ids are the Anax GFF gene ids (via gffread's
  `gene_id`). fate.R's `ag` is `gid()` of the GENESPACE pangenes id; I expect them to match, but
  check the NA count in the join above.
- **Gene IDs vs pangenes**: 05_join.R takes juncea IDs verbatim from `ID=` (`anax_g1`, prefix
  included). They must equal `gid(pangenes$id)` for the filter; strip the prefix in the `sub()` if not.

## Notes

- Both genomes must be **softmasked** (lowercase repeats; checked: Anax 38.4%, Aeshna 47.0%, no N).
  Hard-masked N runs would be read by TOGA2 as assembly gaps. `01_prep_genomes.sh` stops if <1% lowercase.
- Sequence names must not contain dots (`OZ369554.1` -> `OZ369554`). `genomes/*.rename.tsv`
  maps back; the reference GTF is renamed the same way so BED12 and 2bit agree.
- 02 and 04 run everything inside their own allocation (Nextflow local executor, no child jobs).
  `single` partition limits: MaxTime 7 days, MaxNodes 1, MaxMemPerCPU 4000 MB (nodes: 64 CPUs,
  256000 MB; 64 x 4000M is refused). NCPU/NMEM follow `-c` and `--mem-per-cpu` in the 02/04 headers.
- 02 sets `NXF_SYNTAX_PARSER=v1`: Nextflow 26 parses configs strictly by default and rejects the
  `if` block in make_lastz_chains' nextflow.config. TOGA2 (04) runs fine with the default parser.
- 04 writes TOGA2's per-step Nextflow configs to `toga2_nfconfig/` (preprocessing, orthology,
  alignment_3..15); they follow NCPU/NMEM in config.sh.
- SpliceAI is skipped (`--no_spliceai`). It is a GPU/TensorFlow step trained on human splice
  sites; TOGA2 runs without it. To add it: `toga2.py spliceai --help`, then replace
  `--no_spliceai` with `--spliceai_dir <dir>`.
- `--use_long_distance_model` exists for distant pairs; not needed within Aeshnidae.
- The reference annotation goes in as GTF (gffread) because TOGA2's GFF3 parser skips lines
  without `ID=` (BRAKER exon lines) and names each gene after its transcript. `prepare-input`
  derives the isoform table from the GTF and builds reference-specific CESAR2 splice-site
  profiles; the run uses those via `--cesar_profile_dir`.
- Installs: conda env at `/path/to/project/conda/envs/toga` (pip packages from TOGA2's `make` go in it
  too); tools in `/path/to/project/programs/{TOGA2,make_lastz_chains}`; caches (pip, cargo, singularity,
  Nextflow) under `/path/to/project/.cache`, `.singularity`, `.nextflow`.
- Container runtime: SingularityCE 4.4.1 at `/usr/bin/singularity` on login and compute nodes
  (no module); make_lastz_chains runs with `-profile singularity,local`.
- prepare-input wrote 6 of the 8 CESAR2 profiles: none for non-canonical U2 introns (non GT/GC-AG),
  because the Anax gene models contain none. Those two profiles are therefore never applied, but
  `toga2.py run` stops at startup without them, so 04 passes TOGA2's shipped human ones explicitly.
- 04 uses `--annotate_paralogs`: by default TOGA2 aligns a gene's paralogous chains only if it has
  no orthologous chain (`schedulers/preprocess_scheduler.py`), so a neo-Y copy labelled paralogous
  while the neo-X copy is orthologous would get no status. Costs extra alignment jobs genome-wide;
  `--max_chains_per_transcript 10` (default 100) caps them if the run won't fit in 7 days.
- 04 first run: 9.5 h. 82.6% of Anax genes have a functional ortholog (9,755 one-to-one), but 14.8% of
  reference transcripts are called L genome-wide -- probably far more than real losses within one
  family, so read Y calls against `X_call`, not as absolute rates. The run summary's gene-level block
  and its 1:1 sanity warning use denominators that don't match the 15,319 reference genes (39,401 = all
  projections, and 28,507).
- 02 chains (checked with 02b): 95% of Anax CDS aligned; neo-Y 29.0 Mb aligned vs 44.1 Mb on the
  neo-X arm (Y is 167 Mb vs 103 Mb); both m-chromosomes poorly aligned (10-13%).

## 06 -- independent check of the degenerated calls with miniprot (no Anax, no TOGA2)

The A. juncea protein of each gene's neo-X copy is aligned to the neo-Y alone; miniprot reports
frameshifts and in-frame stop codons. Controls in the same run: neo-X proteins of genes whose neo-Y copy
is annotated (positive; sets the gametolog-identity threshold as the 5th percentile of their intact hits)
and ancestral-X proteins (negative; no neo-Y gametolog exists, so their intact high-identity hits are
false positives).

```bash
# one time, login node: add miniprot to the toga env (gffread is already in it)
module load mamba/1.5.8 && mamba install -y -p /path/to/project/conda/envs/toga \
  -c conda-forge -c bioconda --freeze-installed 'miniprot>=0.13'
# laptop: run fate.R (it now writes miniprot_sets.tsv in its working directory), then copy it up
scp miniprot_sets.tsv <user>@<cluster>:/path/to/project/gene_fate/miniprot/
# cluster
cd /path/to/project/gene_fate && sbatch 06_miniprot.slurm      # ~minutes; summary in miniprot.<job>.out
```

Outputs in `miniprot/`: `miniprot_calls.tsv` (per set x Anax gene: best query, class, coverage,
frameshifts, stops, identity, whether the hit lies at a TOGA2 projection of the same gene, and whether it
overlaps an annotated neo-Y gene model), `neoY.miniprot.gff` (raw alignments), `juncea.internal_stop.txt`
(A. juncea models whose own protein has an internal stop; not used as queries).
Classes: intact = >= 90% of the query aligned with no frameshift or stop; broken = >= 50% aligned with
>= 1 frameshift or stop; partial = shorter; none = no neo-Y alignment.
