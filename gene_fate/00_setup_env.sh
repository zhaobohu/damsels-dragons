#!/bin/bash
# Builds TOGA2 and make_lastz_chains into $PROG, skipping what is done (login node: needs internet and the conda env)
set -euo pipefail
source "$(dirname "$(readlink -f "$0")")/config.sh"; activate
command -v gcc >/dev/null || { echo "gcc not on PATH: module load gcc, then rerun"; exit 1; }
mkdir -p "$PROG" "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$GEN" "$CHAINS" "$REFIN"
cd "$PROG"

# 1. TOGA2 (v2.0.10): make compiles, pip-installs into the env, downloads tools and trains models;
#    ld_model.dat is the last file it writes
[ -d TOGA2 ] || git clone --recurse-submodules https://github.com/hillerlab/TOGA2
if [ ! -f TOGA2/models/ld_model.dat ]; then
  (cd TOGA2 && make && python3 toga2.py test -o "$RUN/toga2_test")   # test = local sanity run
fi

# 2. make_lastz_chains (v3.1.7) and its container
[ -d make_lastz_chains ] || git clone https://github.com/hillerlab/make_lastz_chains
[ -f "$MLC/make_lastz_chains.sif" ] || \
  singularity pull "$MLC/make_lastz_chains.sif" docker://ghcr.io/hillerlab/make_lastz_chains:latest

# 3. Nextflow user config: child SLURM jobs get the allocation and the head job's environment
cat > "$NXF_HOME/config" <<CFG
process.clusterOptions = '-A $ALLOC --export=ALL'
CFG
echo "setup done"
