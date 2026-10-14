#!/bin/bash
#SBATCH -D /path/to/project/palaeoptera_synteny/annevo_annotation
#SBATCH -J busco.cpu
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=06:00:00
#SBATCH -a 1-33
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -eox pipefail

mkdir -p logs
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
exec > >(tee -a "logs/busco.${NAME}.out") \
     2> >(grep -v --line-buffered "OpenBLAS Warning" \
          | tee -a "logs/busco.${NAME}.err" >&2)

echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}; ${SLURM_ARRAY_TASK_ID:-0})"
echo "Node: $(hostname)"
echo "Dir: $(pwd)"
echo "Start: $(date)"

export TMPDIR="$(pwd)/tmp.${SLURM_JOB_ID}"
rm -rf "$TMPDIR"
mkdir -p "$TMPDIR"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +x; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
set -x
conda activate /path/to/project/conda/envs/busco

RUN=/path/to/project/palaeoptera_synteny/annevo_annotation
GEN=/path/to/project/data/genomes/${NAME}.fna
GFF=$RUN/gff/${NAME}.gff
AA=$RUN/aa
BUS=$RUN/busco
LIN=$RUN/busco_downloads/lineages/insecta_odb12.2
mkdir -p "$AA" "$BUS"

# skip species already scored (glob, so the odb version is not hardcoded)
compgen -G "$BUS/${NAME}/short_summary.specific.*.txt" >/dev/null \
  && { echo "$NAME already scored"; exit 0; }

# check inputs
[[ -s $GFF ]] || { echo "missing gff: $GFF"; exit 1; }
[[ -s $GEN ]] || { echo "missing genome: $GEN"; exit 1; }
[[ -d $LIN ]] || { echo "missing lineage: $LIN"; exit 1; }
grep -qP '\tCDS\t' "$GFF" || { echo "gff has no CDS: $GFF"; exit 1; }

# longest isoform per gene -> proteins
agat_sp_keep_longest_isoform.pl --gff "$GFF" -o "$AA/${NAME}.longest.gff"
agat_sp_extract_sequences.pl -g "$AA/${NAME}.longest.gff" -f "$GEN" \
  -t cds -p -o "$AA/${NAME}.raw.aa"
sed 's/\*$//' "$AA/${NAME}.raw.aa" > "$AA/${NAME}.aa"
N=$(grep -c '^>' "$AA/${NAME}.aa" || true)
[[ $N -gt 0 ]] || { echo "empty aa: $NAME"; exit 1; }
echo "proteins: $N"

# genes per scaffold, for the record
awk '$3=="gene"{print $1}' "$AA/${NAME}.longest.gff" | sort | uniq -c | sort -k1,1nr

# BUSCO in protein mode, then print the completeness line
busco -i "$AA/${NAME}.aa" -l "$LIN" -m proteins -c 16 \
  --offline --opt-out-run-stats --out_path "$BUS" -o "${NAME}" -f
grep -h -m1 "C:" "$BUS/${NAME}"/short_summary.specific.*.txt
rm -f "$AA/${NAME}.raw.aa"
