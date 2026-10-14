#!/bin/bash
#SBATCH -D /path/to/project/divergence_dating
#SBATCH -J inferences
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Six MCMCtree runs under GTR+G5 (model = 7, chosen by 07.modelsel.sh), then their node ages
set -eox pipefail
shopt -s nullglob

mkdir -p logs
exec > >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.out") \
     2> >(tee -a "logs/${SLURM_JOB_NAME}.${SLURM_JOB_ID}.err" >&2)

echo "Job: ${SLURM_JOB_NAME} (${SLURM_JOB_ID}) on $(hostname)"
echo "Start: $(date)"
rm -rf tmp; mkdir -p tmp
export TMPDIR="$(pwd)/tmp"
trap 'code=$?; echo "Exit: $code at $(date)"' EXIT

{ set +x; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/pcm
set -x

# stop if an input is missing or 07.modelsel.sh has not run
for f in ajun_concat.phy ajun_concat_top75.phy ajun_concat_top50.phy \
         ajun_concat_top25.phy mcmctree.ctl ajun_calib.tre; do
  [ -s "$f" ] || { echo "missing required input: $f" >&2; exit 1; }
done
[ -s modelsel/paml_cp1.iqtree ] || {
  echo "run 07.modelsel.sh first -- it justifies model = 7" >&2; exit 1; }

# two full-data chains, three dS subsets, prior only
M=7
# outdir              alignment              seed model usedata
bash mcmctree.sh full_chain1  ajun_concat.phy        67 $M 2
bash mcmctree.sh full_chain2  ajun_concat.phy        69 $M 2
bash mcmctree.sh ds_top75     ajun_concat_top75.phy  67 $M 2
bash mcmctree.sh ds_top50     ajun_concat_top50.phy  67 $M 2
bash mcmctree.sh ds_top25     ajun_concat_top25.phy  67 $M 2
bash mcmctree.sh prior_only   ajun_concat.phy        67 $M 0

# posterior node ages of every run
printf 'run\tnode\tmean\tet_lo\tet_hi\thpd_lo\thpd_hi\thpd_width\n' \
  > modelsel/ages_all_runs.tsv
for d in full_chain1 full_chain2 ds_top75 ds_top50 ds_top25 prior_only; do
  [ -s "$d/out.txt" ] || continue
  awk -v d="$d" -f ages.awk "$d/out.txt"
done >> modelsel/ages_all_runs.tsv
column -t modelsel/ages_all_runs.tsv
