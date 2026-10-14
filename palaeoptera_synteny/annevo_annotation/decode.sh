#!/bin/bash
#SBATCH -D /path/to/project/palaeoptera_synteny/annevo_annotation
#SBATCH -J decode.cpu
#SBATCH -p single
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --time=08:00:00
#SBATCH -a 1-33
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -eox pipefail

mkdir -p logs
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
exec > >(tee -a "logs/decode.${NAME}.out") \
     2> >(grep -v --line-buffered "OpenBLAS Warning" \
          | tee -a "logs/decode.${NAME}.err" >&2)

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
conda activate /path/to/project/conda/envs/annevo

PRE=/path/to/project/data/genomes/
GEN=${PRE}${NAME}.fna
H5=/path/to/project/palaeoptera_synteny/annevo_annotation/h5/${NAME}.h5
OUT=/path/to/project/palaeoptera_synteny/annevo_annotation/gff
mkdir -p "$OUT"

# skip if already decoded (prot.fa only exists if gffread succeeded on a real gff); check inputs
[[ -s "$OUT/${NAME}.prot.fa" ]] && { echo "$NAME already decoded"; exit 0; }
[[ -s $H5 ]] || { echo "missing h5: $H5"; exit 1; }
[[ -s $GEN ]] || { echo "missing genome: $GEN"; exit 1; }

# decode the ANNEVO predictions into a gff
cd /path/to/project/programs/ANNEVO
python decoding.py -g "$GEN" -p "$H5" -o "$OUT/${NAME}.gff.part" \
  -t 16 --min_intron_length 20 --min_prot_length 100

# decoding.py can exit 0 with a header-only gff (OOM-killed workers), so require CDS records
grep -qP '\tCDS\t' "$OUT/${NAME}.gff.part" || { echo "no CDS decoded: $NAME"; exit 1; }
mv "$OUT/${NAME}.gff.part" "$OUT/${NAME}.gff"

# proteins; the h5 is deleted only after they exist
gffread "$OUT/${NAME}.gff" -g "$GEN" -y "$OUT/${NAME}.prot.fa"
[[ -s "$OUT/${NAME}.prot.fa" ]] || { echo "gffread produced nothing: $NAME"; exit 1; }
rm -f "$H5"
