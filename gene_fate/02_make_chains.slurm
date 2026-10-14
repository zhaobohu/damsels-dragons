#!/bin/bash
#SBATCH -J mlc
#SBATCH -p single
#SBATCH -N 1 -n 1 -c 64 --mem-per-cpu=3500M
#SBATCH -t 7-00:00:00
#SBATCH -o mlc.%j.out
# Whole make_lastz_chains run inside this one job (Nextflow local executor); rerunning resumes
set -euo pipefail
source /path/to/project/gene_fate/config.sh; activate
export NXF_CONTAINER_IMAGE=$MLC/make_lastz_chains.sif
export NXF_SYNTAX_PARSER=v1   # make_lastz_chains' config uses syntax that Nextflow 26's default parser rejects

# job-local resource limits; LASTZ / axtChain / RepeatFiller get 3 GB each instead of the pipeline's
# 24-50 GB (memory and time grow with each retry)
cat > "$CHAINS/local.config" <<CFG
executor {
    cpus   = $NCPU
    memory = ${NMEM}.GB
}
process {
    resourceLimits = [ cpus: $NCPU, memory: ${NMEM}.GB ]
    withName: '.*:LASTZ|.*:REPEAT_FILLER|.*:AXT_CHAIN' {
        cpus   = 1
        memory = { 3.GB * task.attempt }
        time   = { 6.h * task.attempt }
    }
}
CFG

# run the pipeline; chains end up in $CHAINS/results/07_final/
cd "$CHAINS"
nextflow run "$MLC/main.nf" -params-file "$MLC/params.json" \
  -profile singularity,local -c "$CHAINS/local.config" -w "$CHAINS/work" -resume \
  --reference_name "$REF_NAME" --query_name "$QRY_NAME" \
  --reference_genome "$GEN/$REF_NAME.2bit" --query_genome "$GEN/$QRY_NAME.2bit" \
  --outdir "$CHAINS/results"
ls -l "$CHAINS"/results/07_final/
