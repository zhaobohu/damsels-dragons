# Paths and cluster settings, sourced by every script

# edit these
ALLOC=your_slurm_account
NCPU=${SLURM_CPUS_PER_TASK:-64}                             # = '#SBATCH -c' in 02/04
NMEM=$(( NCPU * ${SLURM_MEM_PER_CPU:-3500} / 1024 - 17 ))   # GB for tasks (201 at 64 x 3500M); 17 GB left for Nextflow/TOGA2
RUN=/path/to/project/gene_fate
ENVS=/path/to/project/conda/envs

REF_NAME=Anax_junius            # REFERENCE = annotated genome whose transcripts get projected
QRY_NAME=Aeshna_juncea          # QUERY     = genome being tested (male hap1, GCA_977065675.1)
REF_FA=/path/to/project/data/genomes/Anax.junius.clean.masked.fna         # must be SOFT-masked (lowercase repeats)
QRY_FA=/path/to/project/data/genomes/Aeshna.juncea.clean.masked.fna        # same
REF_GFF=/path/to/project/data/genes/Anax.junius.longest.gff       # reference gene models (mRNA/exon/CDS with ID/Parent)
QRY_GFF=/path/to/project/data/genes/Aeshna.juncea.longest.gff        # only used by 05_join.R

# derived paths, caches and environment
GEN=$RUN/genomes
CHAINS=$RUN/chains
REFIN=$RUN/ref_input
OUT=$RUN/toga2_out
PROG=/path/to/project/programs
TOGA2=$PROG/TOGA2
MLC=$PROG/make_lastz_chains
export SINGULARITY_CACHEDIR=/path/to/project/.singularity/cache   # layer cache for singularity pull
export PIP_CACHE_DIR=/path/to/project/.cache/pip             # TensorFlow etc. during make
export CARGO_HOME=/path/to/project/.cache/cargo              # crates for the Rust build

export NXF_HOME=/path/to/project/.nextflow              # keeps Nextflow out of $HOME quota
export NXF_SINGULARITY_CACHEDIR=/path/to/project/.singularity
export NXF_OPTS="-Xms512m -Xmx4g"   # TOGA2 runs up to 5 Nextflow JVMs at once
export PATH=$TOGA2:$TOGA2/bin:$PATH

activate() {                     # conda env (singularity is at /usr/bin, no module needed)
  local opts=$-; set +u          # conda's shell hooks are not 'set -u' clean
  module load mamba/1.5.8
  eval "$(conda shell.bash hook)"   # scripts don't inherit the conda shell function
  conda activate "$ENVS/toga"
  [[ $opts == *u* ]] && set -u; return 0
}
