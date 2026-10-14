#!/usr/bin/env python3
"""
Per-gene class table (auto / neoX_arm / ancX_arm / neoY / micro, as in the k-mer pipeline) for the
m-chromosome depth check, with shared=1 for genes in 4-way single-copy HOGs (micro and auto genes
also need conserved placement). Writes gene_classes_mchrom.tsv: gene chrom start class shared.

Example:
    python 01_mchrom_gene_sets.py \
      --orthofinder-dir /path/to/project/odonate_synteny/genespace/orthofinder/Results_Jul19 \
      --bed-dir /path/to/project/odonate_synteny/genespace/bed \
      --bed Aeshna_juncea=/path/to/project/gene_fate/depth_check/cds.bed \
      --target-species Aeshna_juncea \
      --other-species Pantala_flavescens Ischnura_elegans Platycnemis_pennipes \
      --micro OZ369553.1 --neox OZ369554.1 --neoy OZ369555.1 \
      --neox-break 103357000 --neox-arm left \
      --autosomes OZ369542.1 OZ369543.1 OZ369544.1 OZ369545.1 OZ369546.1 \
                  OZ369547.1 OZ369548.1 OZ369549.1 OZ369550.1 OZ369551.1 OZ369552.1
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

CLS_ORDER = ("auto", "neoX_arm", "ancX_arm", "neoX", "neoY", "micro")


# BED-like file -> {gene: (chrom, smallest start)}; many rows per gene allowed (cds.bed)
def load_gene_pos(path):
    pos = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            g, c = f[3].strip(), f[0].strip()
            try:
                s = int(f[1])
            except ValueError:
                continue
            if g not in pos or s < pos[g][1]:
                pos[g] = (c, s)
    return pos


# BED for a species: --bed SPECIES=PATH, else <bed-dir>/<species>.bed|.bed.txt|.tsv
def resolve_bed(species, bed_dir, explicit):
    if species in explicit:
        return explicit[species]
    if bed_dir:
        for ext in (".bed", ".bed.txt", ".tsv"):
            p = os.path.join(bed_dir, species + ext)
            if os.path.isfile(p):
                return p
    return None


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orthofinder-dir", required=True)
    p.add_argument("--bed-dir", default=None)
    p.add_argument("--bed", action="append", default=[], metavar="SPECIES=PATH")
    p.add_argument("--target-species", default="Aeshna_juncea")
    p.add_argument("--other-species", nargs="+", required=True)
    p.add_argument("--micro", required=True, help="e.g. OZ369553.1")
    p.add_argument("--neox", default=None, help="e.g. OZ369554.1")
    p.add_argument("--neoy", default=None, help="e.g. OZ369555.1")
    p.add_argument("--neox-break", type=int, default=None,
                   help="fusion breakpoint on the X, bp (e.g. 103357000). "
                        "Omit to treat the X as one class.")
    p.add_argument("--neox-arm", choices=["left", "right"], default="left")
    p.add_argument("--autosomes", nargs="*", default=[],
                   help="Explicit assembled autosomes. Anything not listed and "
                        "not micro/neox/neoy is treated as unplaced and dropped. "
                        "Omit to treat every other sequence as an autosome (not "
                        "recommended if the BED carries unplaced scaffolds).")
    p.add_argument("--exclude-chroms", nargs="*", default=[])
    p.add_argument("--force-partner", action="append", default=[],
                   metavar="SPECIES=CHROM")
    p.add_argument("--strip-re", default=None,
                   help=r"Regex stripped from gene IDs on both sides, e.g. '\.t\d+$'")
    p.add_argument("--out", default="gene_classes_mchrom.tsv")
    p.add_argument("--hog-table", default="mchrom_shared_hogs.tsv")
    args = p.parse_args()

    strip_re = re.compile(args.strip_re) if args.strip_re else None
    norm = (lambda g: strip_re.sub("", g.strip())) if strip_re else (lambda g: g.strip())

    # --bed SPECIES=PATH and --force-partner SPECIES=CHROM pairs
    explicit = {}
    for spec in args.bed:
        if "=" not in spec:
            sys.exit(f"--bed expects SPECIES=PATH, got: {spec}")
        s, pth = spec.split("=", 1)
        explicit[s] = pth
    forced = {}
    for spec in args.force_partner:
        if "=" not in spec:
            sys.exit(f"--force-partner expects SPECIES=CHROM, got: {spec}")
        s, c = spec.split("=", 1)
        forced[s] = c

    target, others = args.target_species, list(args.other_species)
    allsp = [target] + others

    # gene positions: target (chrom, start) and the other species (chrom only)
    tpath = resolve_bed(target, args.bed_dir, explicit)
    if not tpath:
        sys.exit(f"No BED for {target}. Use --bed {target}=/path/to/cds.bed")
    tpos = {norm(k): v for k, v in load_gene_pos(tpath).items()}
    print(f"[bed] {target}: {len(tpos)} genes from {tpath}")

    g2c = {}
    for s in others:
        pth = resolve_bed(s, args.bed_dir, explicit)
        if not pth:
            sys.exit(f"No BED for {s}. Use --bed-dir or --bed {s}=/path")
        g2c[s] = {norm(k): v[0] for k, v in load_gene_pos(pth).items()}
        print(f"[bed] {s}: {len(g2c[s])} genes from {pth}")

    # class of every target gene; excluded and unplaced genes are dropped
    excl = set(args.exclude_chroms)
    autos = set(args.autosomes)
    brk, arm = args.neox_break, args.neox_arm
    named = {args.micro, args.neox, args.neoy} - {None}
    if autos & named:
        sys.exit(f"--autosomes overlaps micro/neox/neoy: {sorted(autos & named)}")

    gclass, dropped = {}, Counter()
    for g, (c, s) in tpos.items():
        if c in excl:
            k = "excluded"
        elif c == args.micro:
            k = "micro"
        elif args.neoy and c == args.neoy:
            k = "neoY"
        elif args.neox and c == args.neox:
            if brk is None:
                k = "neoX"
            else:
                k = "neoX_arm" if ((s < brk) == (arm == "left")) else "ancX_arm"
        elif autos and c not in autos:
            k = "unplaced"
        else:
            k = "auto"
        if k in ("excluded", "unplaced"):
            dropped[k] += 1
        else:
            gclass[g] = k
    print("[class] genes per class: "
          + ", ".join(f"{k}={v}" for k, v in sorted(Counter(gclass.values()).items())))
    if dropped:
        print("[class] dropped: " + ", ".join(f"{k}={v}" for k, v in dropped.items()))
    if not autos:
        print("[class] WARNING: --autosomes not given; every non-named sequence "
              "(including unplaced scaffolds) is in the autosomal baseline")

    # N0.tsv: given directly, or inside the OrthoFinder Results_* directory
    if os.path.isfile(args.orthofinder_dir):
        n0 = args.orthofinder_dir
    else:
        n0 = os.path.join(args.orthofinder_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
        if not os.path.isfile(n0):
            raise FileNotFoundError(f"No N0.tsv found. Looked for:\n  {n0}")
    print(f"[hog] {n0}")

    # HOGs single-copy in all species, with coordinates for every gene
    records = []
    n_total = n_sc = 0
    unmapped = Counter()
    with open(n0, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        missing = [s for s in allsp if s not in header]
        if missing:
            sys.exit("Not in N0.tsv: " + ", ".join(missing)
                     + "\nColumns: " + ", ".join(header[3:]))
        idx = {s: header.index(s) for s in allsp}
        for row in reader:
            n_total += 1
            genes, ok = {}, True
            for s in allsp:
                i = idx[s]
                cell = row[i].strip() if i < len(row) else ""
                gs = [norm(x) for x in cell.split(",") if x.strip()]
                if len(gs) != 1:          # absent OR multi-copy -> not usable
                    ok = False
                    break
                genes[s] = gs[0]
            if not ok:
                continue
            n_sc += 1
            tg = genes[target]
            if tg not in tpos:
                unmapped[target] += 1
                continue
            chrom = {target: tpos[tg][0]}
            bad = False
            for s in others:
                c = g2c[s].get(genes[s])
                if c is None:
                    unmapped[s] += 1
                    bad = True
                    break
                chrom[s] = c
            if bad:
                continue
            records.append({"hog": row[0], "genes": genes, "chrom": chrom})

    print(f"[hog] total {n_total}; single-copy in all {len(allsp)} {n_sc}; "
          f"with coordinates for all {len(records)}")
    for s in allsp:
        if unmapped[s]:
            print(f"[warn] {s}: {unmapped[s]} genes had no BED entry "
                  f"(ID mismatch? try --strip-re)")
    if not records:
        sys.exit("No usable HOGs -- almost certainly a gene-ID namespace mismatch.")

    # modal partner chromosome in each other species, per target chromosome
    partner = defaultdict(lambda: defaultdict(Counter))
    for r in records:
        for s in others:
            partner[r["chrom"][target]][s][r["chrom"][s]] += 1
    modal = {a: {s: partner[a][s].most_common(1)[0][0] for s in others}
             for a in partner}

    # partners of the micro (optionally forced)
    if args.micro not in partner:
        sys.exit(f"No single-copy 4-way HOG places a {target} gene on {args.micro}")
    print(f"\n[partner] modal partner of {args.micro} "
          f"-- check these are clean majorities before trusting the set:")
    for s in others:
        tot = sum(partner[args.micro][s].values())
        for c, n in partner[args.micro][s].most_common(3):
            print(f"    {s}: {c}  {n}/{tot} ({100*n/tot:.0f}%)")
        if s in forced:
            modal[args.micro][s] = forced[s]
            print(f"    {s}: FORCED -> {forced[s]}")

    # shared flag: micro/auto need every ortholog on the modal partner; X arms and neoY do not
    shared, micro_recs = {}, []
    for r in records:
        tg = r["genes"][target]
        k = gclass.get(tg)
        if k is None:
            continue
        if k in ("micro", "auto"):
            a = r["chrom"][target]
            if all(r["chrom"][s] == modal[a][s] for s in others):
                shared[tg] = 1
                if k == "micro":
                    micro_recs.append(r)
        else:
            shared[tg] = 1

    # gene class table, and the shared micro HOGs
    with open(args.out, "w") as fh:
        fh.write("gene\tchrom\tstart\tclass\tshared\n")
        for g, k in sorted(gclass.items()):
            c, s = tpos[g]
            fh.write(f"{g}\t{c}\t{s}\t{k}\t{shared.get(g, 0)}\n")

    with open(args.hog_table, "w") as fh:
        cols = ["HOG", f"{target}_gene", f"{target}_chrom"]
        for s in others:
            cols += [f"{s}_gene", f"{s}_chrom"]
        fh.write("\t".join(cols) + "\n")
        for r in micro_recs:
            out = [r["hog"], r["genes"][target], r["chrom"][target]]
            for s in others:
                out += [r["genes"][s], r["chrom"][s]]
            fh.write("\t".join(out) + "\n")

    # summary
    print("\n[sets] shared=1 genes per class:")
    cnt = Counter(gclass[g] for g in shared)
    tot = Counter(gclass.values())
    for k in CLS_ORDER:
        if k in tot:
            print(f"    {k:9s} {cnt.get(k, 0):6d} shared / {tot[k]:6d} total")
    n_mi = cnt.get("micro", 0)
    if n_mi < 20:
        print(f"\n[warn] only {n_mi} shared micro genes -- the estimate will be noisy; "
              "quote it with the CI, not as a point value")
    print(f"\nwrote {args.out} and {args.hog_table}")


if __name__ == "__main__":
    main()
