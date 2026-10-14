#!/bin/bash
#SBATCH -D /path/to/project/divergence_dating
#SBATCH -J controls
#SBATCH --array=1-500%45
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# Size-matched random controls (09.random_controls.R): task t runs manifest rows t, t+STRIDE, ...,
# NCPU at a time; finished controls are skipped, so it is safe to resubmit
set -eox pipefail
shopt -s nullglob

STRIDE=500            # = array size; keep it when resubmitting a subset of task IDs
FORCE=${FORCE:-0}     # FORCE=1 reruns finished controls too

mkdir -p logs/controls
exec > >(tee -a "logs/controls.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out") \
     2> >(tee -a "logs/controls.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err" >&2)

{ set +x; } 2>/dev/null
module load mamba/1.5.8
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/pcm
set -x

# manifest rows and this task's ID
MAN=controls_manifest.tsv
[ -s "$MAN" ] || { echo "missing $MAN -- run 09.random_controls.R first" >&2; exit 1; }
NROW=$(( $(wc -l < "$MAN") - 1 ))
NCPU=${SLURM_CPUS_PER_TASK:-1}
T=$SLURM_ARRAY_TASK_ID
(( T >= 1 && T <= STRIDE )) || { echo "task $T is outside 1..$STRIDE" >&2; exit 1; }

# complete, and produced after the current alignment was written
done_ok() {
  [ -s "$1/out.txt" ] && grep -q '^t_n15' "$1/out.txt" && [ "$1/out.txt" -nt "$2" ]
}

# run one control, then delete its MCMC trace and PAML scratch (2997 traces would fill the disk)
run_one() {
  local id=$1 phy=$2
  rm -rf "$id"
  bash mcmctree.sh "$id" "$phy" 67 7 2 || return 1
  rm -f "$id/mcmc.txt" "$id/hessian"/tmp*.ctl "$id/hessian"/tmp*.out
  grep -q '^t_n15' "$id/out.txt"
}

# this task's rows
IDS=(); PHYS=()
for (( n = T; n <= NROW; n += STRIDE )); do
  row=$(awk -v n="$n" 'NR==n+1' "$MAN")
  [ -n "$row" ] || { echo "no manifest row $n" >&2; exit 1; }
  id=$(cut -f1 <<< "$row"); phy=$(cut -f3 <<< "$row")
  if [ "$FORCE" != 1 ] && done_ok "$id" "$phy"; then
    echo "[skip] $id already complete"; continue
  fi
  IDS+=("$id"); PHYS+=("$phy")
done
echo "[task $T] ${#IDS[@]} control(s) to run on $NCPU CPU(s): ${IDS[*]:-none}"

# NCPU at a time, one log per control
fail=()
for (( i = 0; i < ${#IDS[@]}; i += NCPU )); do
  pids=()
  for (( j = i; j < i + NCPU && j < ${#IDS[@]}; j++ )); do
    run_one "${IDS[j]}" "${PHYS[j]}" > "logs/controls/${IDS[j]}.log" 2>&1 &
    pids+=($!)
  done
  for (( j = 0; j < ${#pids[@]}; j++ )); do
    wait "${pids[j]}" || fail+=("${IDS[i + j]}")
  done
done

if (( ${#fail[@]} )); then
  echo "[task $T] FAILED: ${fail[*]} -- see logs/controls/<id>.log" >&2
  exit 1
fi
echo "[task $T] done"
