#!/bin/bash
#SBATCH -J repeats
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=16G
#SBATCH -t 01-00:00:00
#SBATCH -D /path/to/project/microchromosome/repeats
#SBATCH -o repeats.%j.out
#SBATCH -e repeats.%j.err
# Repeat landscape (.divsum) and per-hit Kimura table per chromosome from one RepeatMasker .align.
# Expects parse_rm_align.py in the -D directory.
set -euo pipefail

module load mamba/1.5.8
{ set +x; } 2>/dev/null
set +u                      # conda's hook and activate reference unset variables; set -u would abort
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/repeatmodeler
set -u

# inputs
MASKDIR=/path/to/project/annotation/repeats/Aeshna_juncea
ALIGN=$MASKDIR/Aeshna.juncea.clean.fna.align
FAI=$MASKDIR/Aeshna.juncea.clean.fna.fai
UTIL=/path/to/project/conda/envs/repeatmodeler/share/RepeatMasker/util
OUT=.
# every named chromosome, excluding MT (the 16 kb OZ369556.1) and unplaced scaffolds
CHROMS=$(awk '$1 ~ /^OZ/ && $2 >= 1000000 {print $1}' "$FAI")

# check inputs exist
for f in "$ALIGN" "$FAI" "$UTIL/calcDivergenceFromAlign.pl" parse_rm_align.py; do
  [[ -e $f ]] || { echo "ERROR: not found: $f"; exit 1; }
done
mkdir -p $OUT

# sequence sizes, so both landscapes can be expressed as % of sequence
awk '{s += $2} END {print s}' "$FAI"              > $OUT/genome_size.txt
echo "genome $(cat $OUT/genome_size.txt) bp; chromosomes:" $CHROMS

# run both passes in parallel: .divsum summary (no -a: the .align already has Kimura lines) and per-hit table
perl "$UTIL/calcDivergenceFromAlign.pl" -s $OUT/Aeshna.juncea.divsum "$ALIGN" > $OUT/calcdiv.log 2>&1 &
p1=$!
python3 parse_rm_align.py "$ALIGN" $(printf -- "--chrom %s " $CHROMS) 2> $OUT/parse.log \
  | gzip > $OUT/repeats_chroms.tsv.gz &
p2=$!

wait $p1 || { echo "ERROR: calcDivergenceFromAlign.pl failed - see $OUT/calcdiv.log"; exit 1; }
wait $p2 || { echo "ERROR: parse_rm_align.py failed - see $OUT/parse.log"; exit 1; }   # pipefail covers the pipeline

# quick checks of both outputs
echo "--- parser:";      cat $OUT/parse.log
echo "--- summary:";     grep -c "" $OUT/Aeshna.juncea.divsum | xargs echo "lines in divsum:"
grep -q "^Coverage for each repeat class and divergence" $OUT/Aeshna.juncea.divsum \
  || echo "WARNING: divsum has no coverage table - check $OUT/calcdiv.log"

# hits per chromosome, and first rows against the .out (where a header-parsing error would show)
echo "--- hits per chromosome:"
gzip -cd $OUT/repeats_chroms.tsv.gz | awk 'NR > 1 {n[$1]++} END {for (c in n) print c, n[c]}' | sort || true
echo "--- first rows, parser vs .out (start, end and repeat name should agree):"
gzip -cd $OUT/repeats_chroms.tsv.gz | head -4 || true
first=$(echo $CHROMS | awk '{print $1}')
gzip -cd $MASKDIR/Aeshna.juncea.clean.fna.out.gz | awk -v c=$first '$5 == c' | head -3 || true
echo done
