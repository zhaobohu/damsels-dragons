#!/bin/bash
# ModelFinder model choice; exits 1 unless GTR at all three codon positions (later steps use model = 7)
set -eox pipefail
shopt -s nullglob

W=/path/to/project/divergence_dating
cd $W

{ set +x; } 2>/dev/null
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/pcm
set -x

[ -s ajun_concat.phy ]      || { echo "run 06.build_phylip.R first" >&2; exit 1; }
[ -s locus_table_final.tsv ] || { echo "run 05.filter_loci.R first"  >&2; exit 1; }

mkdir -p modelsel
cd modelsel

# split the concatenation into its three codon-position blocks
awk '/^ *[0-9]+ +[0-9]+ *$/ {n++; f=sprintf("cp%d.phy", n)} n {print > f}' \
    ../ajun_concat.phy
for i in 1 2 3; do head -1 cp$i.phy; done

# plain topology for -te: drop the count line and the calibration labels
sed -n '2,$p' ../ajun_calib.tre | sed "s/'[^']*'//g; s/[#@][0-9.]*//g" > topo.tre
cut -c1-120 topo.tre

# ModelFinder per codon position: restricted to the BASEML models with gamma only
# (MCMCtree has no +I or +R), then unrestricted for the supplement
for i in 1 2 3; do
  iqtree3 -s cp$i.phy -te topo.tre -m MF \
          -mset JC,K80,F81,HKY,TN,GTR -mrate E,G \
          --prefix paml_cp$i -T 1 -redo
  iqtree3 -s cp$i.phy -te topo.tre -m MF \
          --prefix free_cp$i -T 1 -redo
done

# per-locus model picks of the final loci, for the supplement
tail -n +2 ../locus_table_final.tsv | cut -f1 | while read -r p; do
  f=../per_locus/$p/$p.iqtree
  [ -s "$f" ] || { echo "missing $f" >&2; continue; }
  printf '%s\t%s\n' "$p" "$(sed -n 's/^Best-fit model[^:]*: //p' "$f" | head -1)"
done > per_locus_models.tsv

# per-locus summary
{
  echo '# per-locus ModelFinder selections'
  cut -f2 per_locus_models.tsv | sed 's/+.*//' | sort | uniq -c | sort -rn |
    awk '{printf "scheme\t%s\t%d\n", $2, $1}'
  printf 'with +I\t\t%d\n' "$(grep -c '+I' per_locus_models.tsv || true)"
  printf 'with +R\t\t%d\n' "$(grep -c '+R' per_locus_models.tsv || true)"
  printf 'HKY85\t\t%d\n'   "$(grep -c 'HKY' per_locus_models.tsv || true)"
} > per_locus_summary.tsv
cat per_locus_summary.tsv

# report, then stop unless GTR was chosen at all three positions
echo '--- restricted to the BASEML model set ---'
sed -n 's/^Best-fit model[^:]*: /  /p' paml_cp1.iqtree paml_cp2.iqtree paml_cp3.iqtree
echo '--- unrestricted ---'
sed -n 's/^Best-fit model[^:]*: /  /p' free_cp1.iqtree free_cp2.iqtree free_cp3.iqtree

n_gtr=$(sed -n 's/^Best-fit model[^:]*: //p' paml_cp1.iqtree paml_cp2.iqtree paml_cp3.iqtree |
        grep -c '^GTR' || true)
if [ "$n_gtr" -ne 3 ]; then
  echo "WARNING: GTR selected at only $n_gtr/3 positions." >&2
  echo "Downstream scripts pass model = 7. Change them before continuing." >&2
  exit 1
fi
echo "GTR selected at all three codon positions -- model = 7 downstream."
