#!/bin/bash
#SBATCH -J miniprot
#SBATCH -p single
#SBATCH -N 1 -n 1 -c 16 --mem-per-cpu=3500M
#SBATCH -t 02:00:00
#SBATCH -o miniprot.%j.out
# miniprot check of the neo-Y "degenerated" calls (neo-X proteins aligned to the neo-Y alone; no Anax, no TOGA2)
# needs $RUN/miniprot/miniprot_sets.tsv from figs/fig4/fig4.R; install miniprot once, on a login node:
#   module load mamba/1.5.8 && mamba install -y -p /path/to/project/conda/envs/toga \
#     -c conda-forge -c bioconda --freeze-installed 'miniprot>=0.13'
set -euo pipefail
source /path/to/project/gene_fate/config.sh
activate
# check tools before anything runs
command -v miniprot >/dev/null || { echo "miniprot not in the toga env: see the install line at the top"; exit 1; }
Rscript -e 'suppressMessages(library(data.table))' ||
  { echo "Rscript/data.table not available: load the R you used for 05_join.R, then resubmit"; exit 1; }
NT=${SLURM_CPUS_PER_TASK:-4}

# work in $RUN/miniprot
MP=$RUN/miniprot; mkdir -p "$MP"; cd "$MP"
Y=OZ369555.1                                   # A. juncea neo-Y, as named in $QRY_FA and $QRY_GFF
[ -s miniprot_sets.tsv ] || { echo "no $MP/miniprot_sets.tsv: copy it from figs/fig4/fig4.R's working directory"; exit 1; }
miniprot --version; gffread --version

# 1. proteins of every A. juncea gene model, one line each, named by gene, terminal stop removed;
#    models with an internal stop are listed, not used
gffread -y juncea.raw.faa -g "$QRY_FA" "$QRY_GFF"
awk '/^>/{if(s!=""){print h; print s}; h=$1; sub(/\.t[0-9]+$/, "", h); s=""; next}
     {s=s $0} END{if(s!=""){print h; print s}}' juncea.raw.faa |
  awk 'NR%2==1{h=$0; next} {sub(/[.*]$/, ""); if ($0 ~ /[.*]/) {print substr(h,2) > "juncea.internal_stop.txt"; next}
       print h; print}' > juncea.faa
touch juncea.internal_stop.txt
echo "proteins: $(grep -c '^>' juncea.faa) usable, $(wc -l < juncea.internal_stop.txt) with internal stops"

# 2. the query proteins (column jg of the sets file)
cut -f4 miniprot_sets.tsv | tail -n +2 | sort -u > query.ids
awk 'NR==FNR{want[$1]; next} /^>/{keep = (substr($1,2) in want)} keep' query.ids juncea.faa > query.faa
echo "queries: $(wc -l < query.ids) requested, $(grep -c '^>' query.faa) found"

# 3. the neo-Y alone as the target, so a neo-X protein cannot simply find its own locus
awk -v id="$Y" '/^>/{p = (substr($1,2) == id)} p' "$QRY_FA" > neoY.fa
[ -s neoY.fa ] || { echo "$Y not found in $QRY_FA"; exit 1; }
echo "target: $Y, $(grep -v '^>' neoY.fa | tr -d '\n' | wc -c) bp"

# 4. align: insect splice model (-j2), up to 200 kb introns (default -G), best alignments only (--outs)
miniprot -t "$NT" -d neoY.mpi neoY.fa
miniprot -t "$NT" -j2 --gff neoY.mpi query.faa > neoY.miniprot.gff
echo "alignments: $(grep -c $'\tmRNA\t' neoY.miniprot.gff)"

# 5. score against the sets and the TOGA2 calls
Rscript "$RUN/06_score_miniprot.R" "$RUN" "$QRY_GFF"
