#!/usr/bin/env python3
"""Species tree from the single-copy orthogroups: MAFFT per locus, concatenation, IQ-TREE."""
import argparse, os, sys, subprocess, csv, glob, re
from collections import defaultdict

# paths - edit here; everything is written under WORKDIR
WORKDIR = "/path/to/project/palaeoptera_synteny/pairwise_synteny"
OFD     = "/path/to/project/palaeoptera_synteny/genespace_run/genespace/orthofinder/Results_Aug18"
OUTDIR  = "tree"          # relative to WORKDIR
MAFFT   = "mafft"
IQTREE  = "iqtree2"

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=WORKDIR)
ap.add_argument("--ofd", default=OFD, help="OrthoFinder Results_* directory")
ap.add_argument("--out", default=OUTDIR, help="output directory, relative to --workdir")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--resume", action="store_true", help="omit -redo so IQ-TREE restarts from its checkpoint; use for a long sbatch job that may hit walltime")
ap.add_argument("--bootstrap", type=int, default=1000, help="ultrafast bootstraps; 0 to skip")
ap.add_argument("--model", default="MFP", help="IQ-TREE model; MFP = ModelFinder, or e.g. LG+F+G4 for a quick fixed-model tree")
ap.add_argument("--min-taxa", type=int, default=29, help="drop orthogroups with fewer sequences")
ap.add_argument("--iqtree", default=IQTREE)
ap.add_argument("--mafft", default=MAFFT)
a = ap.parse_args()

# resolve the OrthoFinder path before chdir; make the output directories
a.ofd = os.path.abspath(a.ofd)
os.makedirs(a.workdir, exist_ok=True)
os.chdir(a.workdir)

ogd  = os.path.join(a.ofd, "Orthogroups")
scos = os.path.join(a.ofd, "Single_Copy_Orthologue_Sequences")
for p in (ogd, scos):
    if not os.path.isdir(p):
        sys.exit(f"missing directory: {p}")
os.makedirs(os.path.join(a.out, "aln"), exist_ok=True)
os.makedirs(os.path.join(a.out, "iqtree"), exist_ok=True)

# 1. gene -> species, from Orthogroups.tsv (OrthoFinder -X: headers carry no species)
gene2sp = {}
with open(os.path.join(ogd, "Orthogroups.tsv")) as fh:
    rd = csv.reader(fh, delimiter="\t")
    header = next(rd); species = header[1:]
    for row in rd:
        for sp, cell in zip(species, row[1:]):
            for g in cell.split(","):
                g = g.strip()
                if g:
                    gene2sp[g] = sp
print(f"[map] {len(gene2sp)} genes -> {len(species)} species", file=sys.stderr)


def read_fasta(f):
    name, seq, out = None, [], []
    for line in open(f):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def aln_ok(path, n_expected):
    """a usable alignment: parses, has the expected sequences, equal lengths, > 0 sites"""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    try:
        recs = read_fasta(path)
    except Exception:
        return False
    if len(recs) != n_expected:
        return False
    L = {len(seq) for _, seq in recs}
    return len(L) == 1 and L.pop() > 0


# 2. rename headers to species and align each orthogroup with MAFFT
# (a cached alignment is reused only if it passes aln_ok; a locus MAFFT cannot align is dropped)
files = sorted(glob.glob(os.path.join(scos, "*.fa")) + glob.glob(os.path.join(scos, "*.fasta")))
if not files:
    sys.exit(f"no fasta files in {scos}")
kept, dropped, failed = [], 0, 0
for f in files:
    og = re.sub(r"\.(fa|fasta)$", "", os.path.basename(f))
    recs = read_fasta(f)
    named = []
    for gname, seq in recs:
        sp = gene2sp.get(gname)
        if sp is None:
            print(f"[warn] {og}: gene {gname} not in Orthogroups.tsv; skipped", file=sys.stderr); continue
        named.append((sp, seq))

    # drop loci with too few taxa, duplicate species or empty sequences
    if len(named) < a.min_taxa or len({s for s, _ in named}) != len(named):
        dropped += 1; continue
    named = [(sp, seq) for sp, seq in named if seq.strip()]
    if len(named) < a.min_taxa:
        dropped += 1; continue

    tmp = os.path.join(a.out, "aln", og + ".in.fa")
    aln = os.path.join(a.out, "aln", og + ".aln.fa")
    if not aln_ok(aln, len(named)):
        if os.path.exists(aln):
            print(f"[redo] {og}: cached alignment unusable, realigning", file=sys.stderr)
            os.remove(aln)
        with open(tmp, "w") as fh:
            for sp, seq in named:
                fh.write(f">{sp}\n{seq}\n")
        part = aln + ".part"
        cmd = [a.mafft, "--quiet", "--thread", str(a.threads)]
        cmd += ["--auto"]
        cmd.append(tmp)
        with open(part, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
        os.remove(tmp)
        if not aln_ok(part, len(named)):
            os.remove(part)
            print(f"[warn] {og}: MAFFT produced an unusable alignment; locus dropped", file=sys.stderr)
            failed += 1; continue
        os.replace(part, aln)          # only a validated file gets the final name
    kept.append((og, aln))
print(f"[align] {len(kept)} orthogroups aligned, {dropped} dropped (< {a.min_taxa} taxa, "
      f"duplicate species or empty sequences), {failed} failed alignment", file=sys.stderr)
if not kept:
    sys.exit("nothing to concatenate")

# 3. concatenate into a supermatrix (missing species as gaps) and write a partition file
concat = defaultdict(list); pos = 1; parts = []
bad = 0
for og, aln in kept:
    recs = dict(read_fasta(aln))
    if not recs:
        print(f"[warn] {og}: empty alignment at concatenation; skipped", file=sys.stderr)
        bad += 1; continue
    L = {len(v) for v in recs.values()}
    if len(L) != 1:
        print(f"[warn] {og}: ragged alignment ({sorted(L)}); skipped", file=sys.stderr)
        bad += 1; continue
    L = L.pop()
    for sp in species:
        concat[sp].append(recs.get(sp, "-" * L))
    parts.append(f"LG, {og} = {pos}-{pos + L - 1}"); pos += L
if bad:
    print(f"[warn] {bad} locus/loci skipped at concatenation; delete tree/aln and rerun "
          "if that number is large", file=sys.stderr)
if not parts:
    sys.exit("no usable alignments; delete the aln/ directory and rerun")
cfa = os.path.join(a.out, "concat.faa")
with open(cfa, "w") as fh:
    for sp in species:
        fh.write(f">{sp}\n{''.join(concat[sp])}\n")
with open(os.path.join(a.out, "partitions.txt"), "w") as fh:
    fh.write("\n".join(parts) + "\n")
print(f"[concat] {len(species)} taxa x {pos - 1} sites, {len(parts)} partitions -> {cfa}", file=sys.stderr)

# 4. IQ-TREE on the whole supermatrix (one model; partitions.txt is not used)
pre = os.path.join(a.out, "iqtree", "tree")
cmd = [a.iqtree, "-s", cfa, "-T", str(a.threads), "--prefix", pre]
if not a.resume:
    cmd.append("-redo")   # -redo discards any checkpoint and starts clean
if a.bootstrap > 0:
    cmd += ["-B", str(a.bootstrap)]
cmd += ["-m", a.model]
print("[iqtree] " + " ".join(cmd), file=sys.stderr)
subprocess.run(cmd, check=True)
print(f"[done] {pre}.treefile", file=sys.stderr)
