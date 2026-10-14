#!/bin/bash
#SBATCH -D /path/to/project/palaeoptera_synteny/annevo_annotation
#SBATCH -J predict.gpu
#SBATCH -p gpu2
#SBATCH --gres=gpu:1
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time=1-00:00:00
#SBATCH -a 1-33%3
#SBATCH -o /dev/null
#SBATCH -e /dev/null

set -eox pipefail

mkdir -p logs
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" species.txt)
exec > >(tee -a "logs/predict.${NAME}.out") \
     2> >(grep -v --line-buffered "OpenBLAS Warning" \
          | tee -a "logs/predict.${NAME}.err" >&2)

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
# one-time setup, on the login node:
#   git clone https://github.com/xjtu-omics/ANNEVO
#   conda create -y -p /path/to/project/conda/envs/annevo python=3.10
#   conda activate /path/to/project/conda/envs/annevo
#   cd /path/to/project/programs/ANNEVO && pip install -r requirements.txt
#   (download the Invertebrate model to $MOD per the repo README)

# check that a GPU is visible
python -c "import torch; assert torch.cuda.is_available(), 'no CUDA'; print(torch.cuda.get_device_name(0))"

PRE=/path/to/project/data/genomes/
GEN=${PRE}${NAME}.fna
mkdir -p /path/to/project/palaeoptera_synteny/annevo_annotation/h5
H5=/path/to/project/palaeoptera_synteny/annevo_annotation/h5/${NAME}.h5
MOD=/path/to/project/programs/ANNEVO/saved_model/ANNEVO_Insecta.pt

# skip if already predicted; check the genome
[[ -s $H5 ]] && { echo "$NAME already predicted"; exit 0; }
[[ -s $GEN ]] || { echo "missing genome: $GEN"; exit 1; }

# ANNEVO prediction on the GPU -> h5 (decoded into a gff by decode.sh)
cd /path/to/project/programs/ANNEVO
python prediction.py -g "$GEN" -m "$MOD" -l Insecta -s 200 \
  -p "${H5}.part" --batch_size 16 --num_workers 8
mv "${H5}.part" "$H5"
