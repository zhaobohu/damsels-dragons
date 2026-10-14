#!/bin/bash
# Per-dataset Hessian, then MCMCtree; called by 08.inferences.sbatch and 10.controls.sbatch
# usage: bash mcmctree.sh <outdir> <alignment.phy> <seed> [model, 7 = GTR] [usedata, 2 or 0 = prior]
set -eox pipefail

{ set +x; } 2>/dev/null
eval "$(conda shell.bash hook)"
conda activate /path/to/project/conda/envs/pcm
set -x

W=/path/to/project/divergence_dating
O=$1; PHY=$2; SEED=${3:--1}; MODEL=${4:-}; USE=${5:-2}

# copy the template ctl, tree and alignment into <outdir> and set seqfile, seed, model
cd $W
[ -s "$PHY" ] || { echo "missing $PHY" >&2; exit 1; }
[ -e "$O" ] && { echo "$O exists, refusing to overwrite" >&2; exit 1; }

mkdir -p "$O"
cp mcmctree.ctl ajun_calib.tre "$O/"
cp "$PHY" "$O/"
sed -i "s|^ *seqfile *=.*|seqfile = $(basename "$PHY")|" "$O/mcmctree.ctl"
sed -i "s|^ *seed *=.*|seed = $SEED|" "$O/mcmctree.ctl"
[ -n "$MODEL" ] && sed -i "s|^ *model *=.*|model = $MODEL|" "$O/mcmctree.ctl"

# prior only; otherwise the Hessian (usedata = 3) in hessian/, then the approximate likelihood
if [ "$USE" = 0 ]; then
  sed -i 's/^ *usedata *=.*/usedata = 0/' "$O/mcmctree.ctl"
  cd "$O"
  mcmctree mcmctree.ctl
else
  mkdir -p "$O/hessian"
  cp "$O/mcmctree.ctl" "$O/ajun_calib.tre" "$O/$(basename "$PHY")" "$O/hessian/"
  sed -i 's/^ *usedata *=.*/usedata = 3/' "$O/hessian/mcmctree.ctl"

  cd "$O/hessian"
  mcmctree mcmctree.ctl
  [ -s out.BV ] || { echo "no out.BV produced" >&2; exit 1; }
  mv out.BV ../in.BV

  cd ..
  # usedata = 2 needs the file name, or MCMCtree aborts with "error: file name empty."
  sed -i 's/^ *usedata *=.*/usedata = 2 in.BV/' mcmctree.ctl
  mcmctree mcmctree.ctl
fi

# md5s of the Hessian and alignment, and the ctl settings used
: > provenance.md5
[ -s in.BV ] && md5sum in.BV >> provenance.md5
md5sum "$(basename "$PHY")" >> provenance.md5
grep -E '^ *(seqfile|usedata|model|alpha|ncatG|clock|seed)' mcmctree.ctl >> provenance.md5
