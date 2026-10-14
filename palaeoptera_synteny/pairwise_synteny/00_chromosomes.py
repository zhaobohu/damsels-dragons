#!/usr/bin/env python3
"""Chromosome set per genome from NCBI sequence reports -> assemblies.tsv, chromosomes.tsv."""
import argparse, json, os, subprocess, sys, glob, collections

# paths - edit here
WORKDIR   = "/path/to/project/palaeoptera_synteny/pairwise_synteny"
BED_DIR   = "/path/to/project/palaeoptera_synteny/genespace_run/genespace/bed"
DOWNLOAD  = "/path/to/project/palaeoptera_synteny/download_genomes"          # from download.sbatch
DATASETS  = "/path/to/project/palaeoptera_synteny/download_genomes/datasets"  # the CLI binary
DATA_REPORT = os.path.join(DOWNLOAD, "raw/ncbi_dataset/data/assembly_data_report.jsonl")
SEQREP_DIR  = "seqrep"          # relative to WORKDIR
OUT_ASM     = "assemblies.tsv"  # relative to WORKDIR
OUT_CHROM   = "chromosomes.tsv" # relative to WORKDIR
API_KEY     = os.environ.get("NCBI_API_KEY", "")

ap = argparse.ArgumentParser()
ap.add_argument("--workdir", default=WORKDIR)
ap.add_argument("--bed", default=BED_DIR)
ap.add_argument("--data-report", default=DATA_REPORT)
ap.add_argument("--datasets", default=DATASETS)
ap.add_argument("--seqrep", default=SEQREP_DIR)
ap.add_argument("--skip-download", action="store_true")
a = ap.parse_args()

# resolve inputs before chdir, so relative arguments keep working
a.bed = os.path.abspath(a.bed)
a.data_report = os.path.abspath(a.data_report)
if os.path.exists(a.datasets) or os.sep in a.datasets:
    a.datasets = os.path.abspath(a.datasets)
os.makedirs(a.workdir, exist_ok=True)
os.chdir(a.workdir)
seqrep = a.seqrep
os.makedirs(seqrep, exist_ok=True)

# 1. genomes = bed basenames
beds = sorted(os.path.basename(f)[:-4] for f in glob.glob(os.path.join(a.bed, "*.bed")))
if not beds:
    sys.exit(f"no .bed files in {a.bed}")
print(f"[bed] {len(beds)} genomes", file=sys.stderr)

if not os.path.exists(a.data_report):
    sys.exit(f"missing {a.data_report}\n"
             "  that file comes from the `datasets download genome ... --dehydrated` step;\n"
             "  if the download directory was cleaned, regenerate it with\n"
             "  `datasets summary genome taxon 33339 --assembly-level chromosome,complete "
             "--assembly-source genbank --as-json-lines > assembly_data_report.jsonl`")

# (organism, accession) for every assembly in the download's data report
pairs = []
with open(a.data_report) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        acc = r.get("accession") or r.get("currentAccession")
        o = r.get("organism") or {}
        name = (o.get("organismName") or r.get("organismName") or "").replace(" ", "_")
        if acc and name:
            pairs.append((name, acc))

by_name = collections.defaultdict(list)
for name, acc in pairs:
    by_name[name].append(acc)

# candidate accessions per genome (several = haplotype pairs; chosen by agreement with the bed in step 3)
cands, missing = {}, []
for g in beds:
    c = sorted(set(by_name.get(g, [])))
    if not c:
        # genus_species prefix, in case the bed name carries a strain suffix
        c = sorted({acc for n, acc in pairs if g.startswith(n) or n.startswith(g)})
    if not c:
        missing.append(g)
    else:
        cands[g] = c
if missing:
    print(f"[warn] no accession in the data report for: {', '.join(missing)}", file=sys.stderr)
multi = {g: c for g, c in cands.items() if len(c) > 1}
if multi:
    print(f"[info] {len(multi)} genome(s) have several assemblies (haplotype pairs); "
          "choosing by agreement with the bed", file=sys.stderr)


def seqrep_path(acc):
    hits = glob.glob(os.path.join(seqrep, "**", acc, "sequence_report.jsonl"), recursive=True)
    return hits[0] if hits else None


# 2. fetch the missing sequence reports (metadata only) with the datasets CLI
all_acc = sorted({acc for c in cands.values() for acc in c})
need = [acc for acc in all_acc if seqrep_path(acc) is None]
if need and not a.skip_download:
    listfile = os.path.join(seqrep, "accessions.txt")
    with open(listfile, "w") as fh:
        fh.write("\n".join(need) + "\n")
    zipf = os.path.join(seqrep, "seqrep.zip")
    cmd = [a.datasets, "download", "genome", "accession", "--inputfile", listfile,
           "--include", "seq-report", "--no-progressbar", "--filename", zipf]
    if API_KEY:
        cmd += ["--api-key", API_KEY]
    print("[datasets] " + " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True)
    subprocess.run(["unzip", "-o", "-q", zipf, "-d", seqrep], check=True)
elif need:
    sys.exit(f"{len(need)} sequence reports missing and --skip-download was given")


def get(rec, *keys):
    for k in keys:
        if k in rec and rec[k] not in (None, ""):
            return rec[k]
    return None


# 3. per genome, keep the candidate assembly with the most chromosomes in the bed
asm, rows, report = [], [], []
for g in beds:
    if g not in cands:
        continue

    # genes per scaffold in the bed
    counts = collections.Counter()
    with open(os.path.join(a.bed, g + ".bed")) as fh:
        for line in fh:
            counts[line.split("\t")[0]] += 1

    # chromosomes (assembled-molecule + Chromosome) of each candidate that appear in the bed
    scored = []
    for acc in cands[g]:
        p = seqrep_path(acc)
        if p is None:
            print(f"[warn] no sequence report for {g} ({acc})", file=sys.stderr)
            continue
        ch, roles = [], collections.Counter()
        with open(p) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                role = get(r, "role", "sequenceRole")
                mtype = get(r, "assignedMoleculeLocationType", "assignedMoleculeLocationOrType")
                if role is None:
                    sys.exit(f"{p}: no role field; keys present: {sorted(r)}")
                roles[(role, mtype)] += 1
                if role != "assembled-molecule" or mtype != "Chromosome":
                    continue
                gb = get(r, "genbankAccession"); rs = get(r, "refseqAccession")
                chrom = gb if gb in counts else (rs if rs in counts else None)
                if chrom is None:
                    continue
                ch.append((chrom, get(r, "chrName", "sequenceName") or "", counts[chrom]))
        scored.append((len(ch), acc, ch, roles))
    if not scored:
        print(f"[warn] {g}: no sequence reports at all - skipped", file=sys.stderr)
        continue

    scored.sort(key=lambda x: -x[0])
    kept, acc, ch, roles = scored[0]
    if kept == 0:
        print(f"[warn] {g}: none of {'/'.join(cands[g])} has chromosomes matching the bed "
              "- skipped. Check that the bed was built from one of these assemblies.",
              file=sys.stderr)
        continue
    if len(scored) > 1:
        alt = ", ".join(f"{a}:{n}" for n, a, _, _ in scored[1:])
        print(f"[choose] {g}: {acc} ({kept} chromosomes match the bed; other candidates {alt})",
              file=sys.stderr)
    asm.append((g, acc))
    for chrom, name, n in ch:
        rows.append((g, chrom, name, n))
    report.append((g, kept, roles))

# write the chosen assemblies and the chromosome table
with open(OUT_ASM, "w") as fh:
    fh.write("genome\taccession\n")
    for g, acc in asm:
        fh.write(f"{g}\t{acc}\n")
print(f"[wrote] {OUT_ASM}  ({len(asm)} of {len(beds)} genomes)", file=sys.stderr)
if len(asm) != len(beds):
    print("[warn] some genomes are missing from assemblies.tsv; fix before step 2, "
          "which requires bed/ and chromosomes.tsv to cover the same genomes", file=sys.stderr)

if not rows:
    sys.exit("no chromosomes selected; check the sequence report fields")

with open(OUT_CHROM, "w") as fh:
    fh.write("genome\tchr\tmolecule\tgenes\n")
    for r in sorted(rows, key=lambda x: (x[0], -x[3])):
        fh.write("\t".join(str(v) for v in r) + "\n")

# summary per genome: chromosomes kept, their gene counts, and the sequence roles seen
print(f"\n[wrote] {OUT_CHROM}  ({len(rows)} chromosomes across {len(report)} genomes)\n",
      file=sys.stderr)
print(f"{'genome':28s} {'chr':>4s} {'min_genes':>10s} {'max_genes':>10s}  roles seen", file=sys.stderr)
for g, kept, roles in report:
    gg = [r[3] for r in rows if r[0] == g]
    mn = f"{min(gg):10d}" if gg else f"{'-':>10s}"
    mx = f"{max(gg):10d}" if gg else f"{'-':>10s}"
    other = ", ".join(f"{k[0]}/{k[1]}:{v}" for k, v in sorted(roles.items()))
    print(f"{g:28s} {kept:4d} {mn} {mx}  {other}", file=sys.stderr)
