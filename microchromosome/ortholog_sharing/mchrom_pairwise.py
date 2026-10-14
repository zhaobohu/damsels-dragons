#!/usr/bin/env python3
"""
Pairwise m-chromosome homology check: for each pair of species in --micro, HOGs with >=1 gene on
both m-chromosomes (observed) vs marginal_A * marginal_B / total_HOGs (expected if independent).

Example:
  python mchrom_pairwise.py \\
    --orthofinder-dir <PALAEO_RESULTS> --bed-dir <PALAEO_BED_DIR> \\
    --bed Aeshna_juncea=../cds_depth/cds.annevo.bed --strip-re '\\.t\\d+$' \\
    --micro Aeshna_juncea=OZ369553.1 --micro Ischnura_elegans=NC_060258.1 \\
    --micro Platycnemis_pennipes=OW121858.1 --micro Tanypteryx_hageni=CM054119.1
"""
import argparse
import csv
import itertools
import os
import re
import sys
from collections import Counter


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orthofinder-dir", required=True)
    p.add_argument("--bed-dir", default=None)
    p.add_argument("--bed", action="append", default=[], metavar="SPECIES=PATH")
    p.add_argument("--micro", action="append", required=True, metavar="SPECIES=CHROM")
    p.add_argument("--strip-re", default=None)
    args = p.parse_args()

    strip_re = re.compile(args.strip_re) if args.strip_re else None
    norm = (lambda g: strip_re.sub("", g.strip())) if strip_re else (lambda g: g.strip())
    explicit = dict(s.split("=", 1) for s in args.bed)
    micro = dict(s.split("=", 1) for s in args.micro)
    species = list(micro)

    # gene -> chromosome per species, from --bed SPECIES=PATH or <bed-dir>/<species>.bed|.bed.txt|.tsv
    g2c = {}
    for s in species:
        pth = None
        if s in explicit:
            pth = explicit[s]
        elif args.bed_dir:
            for ext in (".bed", ".bed.txt", ".tsv"):
                cand = os.path.join(args.bed_dir, s + ext)
                if os.path.isfile(cand):
                    pth = cand
                    break
        if not pth:
            sys.exit(f"No BED for {s}")
        m = {}
        with open(pth) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    m[f[3].strip()] = f[0].strip()
        g2c[s] = {norm(k): v for k, v in m.items()}
        n_on = sum(1 for c in g2c[s].values() if c == micro[s])
        print(f"[bed] {s}: {len(g2c[s])} genes, {n_on} on {micro[s]}")
        if n_on == 0:
            sys.exit(f"  {micro[s]} carries no genes in {s}'s BED")

    # N0.tsv: given directly, or inside the OrthoFinder Results_* directory
    if os.path.isfile(args.orthofinder_dir):
        n0 = args.orthofinder_dir
    else:
        n0 = os.path.join(args.orthofinder_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
        if not os.path.isfile(n0):
            raise FileNotFoundError(f"No N0.tsv found. Looked for:\n  {n0}")

    # per species: HOGs with >=1 gene on its m-chromosome, and HOGs with >=1 gene anywhere
    total = 0
    on_m = {s: set() for s in species}          # HOG ids with >=1 gene on the m
    any_present = {s: set() for s in species}   # HOG ids with >=1 gene anywhere
    with open(n0, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        missing = [s for s in species if s not in header]
        if missing:
            sys.exit("Not in N0.tsv: " + ", ".join(missing))
        idx = {s: header.index(s) for s in species}
        for row in reader:
            total += 1
            hid = row[0]
            for s in species:
                i = idx[s]
                cell = row[i].strip() if i < len(row) else ""
                genes = [norm(x) for x in cell.split(",") if x.strip()]
                if genes:
                    any_present[s].add(hid)
                if any(g2c[s].get(g) == micro[s] for g in genes):
                    on_m[s].add(hid)

    # marginal counts
    print(f"\n[hog] total HOGs {total}")
    print("\n[marginal] HOGs with >=1 gene on the m-chromosome:")
    for s in species:
        print(f"    {s:24s} {len(on_m[s]):6d}   (present anywhere: {len(any_present[s])})")

    # pairwise observed vs expected
    print("\n[pairwise] observed / expected-if-independent / enrichment:")
    print(f"    {'':24s}" + "".join(f"{s[:14]:>18s}" for s in species))
    for a in species:
        cells = []
        for b in species:
            if a == b:
                cells.append(f"{'—':>18s}")
                continue
            obs = len(on_m[a] & on_m[b])
            exp = len(on_m[a]) * len(on_m[b]) / total if total else float("nan")
            enr = obs / exp if exp else float("inf")
            cells.append(f"{obs:5d}/{exp:6.1f} ={enr:5.1f}x")
        print(f"    {a:24s}" + "".join(cells))

    # of the HOGs on A's m that B has at all, the fraction on B's m
    print("\n[conditional] of HOGs on A's m-chromosome that B has at all, fraction on B's m-chromosome:")
    print(f"    {'A \\ B':24s}" + "".join(f"{s[:14]:>18s}" for s in species))
    for a in species:
        cells = []
        for b in species:
            if a == b:
                cells.append(f"{'—':>18s}")
                continue
            shared_any = on_m[a] & any_present[b]
            k = len(shared_any & on_m[b])
            cells.append(f"{k:5d}/{len(shared_any):5d} ={100*k/len(shared_any) if shared_any else 0:5.1f}%")
        print(f"    {a:24s}" + "".join(cells))

    # HOGs on the m in every species, and in every subset of species
    if len(species) >= 3:
        allm = set.intersection(*(on_m[s] for s in species))
        print(f"\n[all] HOGs on the m-chromosome in every listed species: {len(allm)}")
        for k in range(len(species) - 1, 1, -1):
            for combo in itertools.combinations(species, k):
                inter = set.intersection(*(on_m[s] for s in combo))
                print(f"    {' + '.join(c.split('_')[0] for c in combo):45s} {len(inter):6d}")


if __name__ == "__main__":
    main()
