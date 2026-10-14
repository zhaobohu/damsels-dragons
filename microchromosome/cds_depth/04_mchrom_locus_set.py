#!/usr/bin/env python3
"""
m-chromosome locus set from one OrthoFinder run, with the m-chromosome named for every species.
  denominator  HOGs in which every comparison species has >=1 gene on its m-chromosome
  numerator    of those, HOGs with exactly one target gene, on the target m-chromosome
Also writes gene_classes_mchrom.tsv (gene chrom start class shared) for 03_mchrom_ratio.R.

Example:
  python 04_mchrom_locus_set.py \
    --orthofinder-dir /path/to/palaeoptera/orthofinder/Results_XXX \
    --bed-dir /path/to/palaeoptera/genespace/bed \
    --target-species Aeshna_juncea \
    --micro Aeshna_juncea=OZ369553.1 \
    --micro Ischnura_elegans=NC_060258.1 \
    --micro Platycnemis_pennipes=OW121858.1 \
    --micro Tanypteryx_hageni=CM054119.1 \
    --neox OZ369554.1 --neoy OZ369555.1 --neox-break 103356500 --neox-arm left \
    --autosomes OZ369542.1 ... OZ369552.1
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

CLS_ORDER = ("auto", "neoX_arm", "ancX_arm", "neoX", "neoY", "micro")


# SPECIES=VALUE arguments -> dict
def kv_list(items, flag):
    out = {}
    for spec in items:
        if "=" not in spec:
            sys.exit(f"{flag} expects SPECIES=VALUE, got: {spec}")
        k, v = spec.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--orthofinder-dir", required=True)
    p.add_argument("--bed-dir", default=None)
    p.add_argument("--bed", action="append", default=[], metavar="SPECIES=PATH")
    p.add_argument("--target-species", required=True)
    p.add_argument("--micro", action="append", required=True, metavar="SPECIES=CHROM",
                   help="m-chromosome ID per species; repeat for target + each comparison species")
    p.add_argument("--neox", default=None)
    p.add_argument("--neoy", default=None)
    p.add_argument("--neox-break", type=int, default=None)
    p.add_argument("--neox-arm", choices=["left", "right"], default="left")
    p.add_argument("--autosomes", nargs="*", default=[])
    p.add_argument("--exclude-chroms", nargs="*", default=[])
    p.add_argument("--min-comp", type=int, default=None,
                   help="minimum number of comparison species that must have >=1 gene on "
                        "their m-chromosome for a HOG to enter the denominator "
                        "(default: all of them)")
    p.add_argument("--strip-re", default=None)
    p.add_argument("--out", default="gene_classes_mchrom.tsv")
    p.add_argument("--hog-table", default="mchrom_locus_hogs.tsv")
    args = p.parse_args()

    strip_re = re.compile(args.strip_re) if args.strip_re else None
    norm = (lambda g: strip_re.sub("", g.strip())) if strip_re else (lambda g: g.strip())

    explicit = kv_list(args.bed, "--bed")
    micro = kv_list(args.micro, "--micro")
    target = args.target_species
    if target not in micro:
        sys.exit(f"--micro must include the target species ({target})")
    comps = [s for s in micro if s != target]
    if not comps:
        sys.exit("need at least one comparison species in --micro")
    allsp = [target] + comps
    tmicro = micro[target]
    min_comp = len(comps) if args.min_comp is None else args.min_comp
    if not 1 <= min_comp <= len(comps):
        sys.exit(f"--min-comp must be between 1 and {len(comps)}")

    # gene -> (chrom, smallest start) per species; BED from --bed or <bed-dir>/<species>.bed|.bed.txt|.tsv
    # stop if a species' m-chromosome carries no genes (ID namespace check)
    pos = {}
    for s in allsp:
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
            sys.exit(f"No BED for {s}. Use --bed-dir or --bed {s}=/path")
        raw = {}
        with open(pth) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 4:
                    continue
                g, c = f[3].strip(), f[0].strip()
                try:
                    start = int(f[1])
                except ValueError:
                    continue
                if g not in raw or start < raw[g][1]:
                    raw[g] = (c, start)
        pos[s] = {norm(k): v for k, v in raw.items()}
        chroms = Counter(v[0] for v in pos[s].values())
        n_on = chroms.get(micro[s], 0)
        print(f"[bed] {s}: {len(pos[s])} genes; {n_on} on m-chromosome {micro[s]}  ({pth})")
        if n_on == 0:
            top = ", ".join(f"{c}({n})" for c, n in chroms.most_common(8))
            sys.exit(f"  {micro[s]} carries no genes in {s}'s BED. Chromosome IDs seen: {top} ...\n"
                     "  Check the accession namespace (RefSeq NC_/NW_ vs GenBank CM_/OW_) for this run.")

    # class of every target gene; excluded and unplaced genes are dropped
    excl, autos = set(args.exclude_chroms), set(args.autosomes)
    brk, arm = args.neox_break, args.neox_arm
    named = {tmicro, args.neox, args.neoy} - {None}
    if autos & named:
        sys.exit(f"--autosomes overlaps micro/neox/neoy: {sorted(autos & named)}")

    gclass, dropped = {}, Counter()
    for g, (c, s) in pos[target].items():
        if c in excl:
            k = "excluded"
        elif c == tmicro:
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
    print("[class] target genes per class: "
          + ", ".join(f"{k}={v}" for k, v in sorted(Counter(gclass.values()).items())))
    if dropped:
        print("[class] dropped: " + ", ".join(f"{k}={v}" for k, v in dropped.items()))
    if not autos:
        print("[class] WARNING: --autosomes not given; unplaced scaffolds are in the autosomal baseline")

    # N0.tsv: given directly, or inside the OrthoFinder Results_* directory
    if os.path.isfile(args.orthofinder_dir):
        n0 = args.orthofinder_dir
    else:
        n0 = os.path.join(args.orthofinder_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
        if not os.path.isfile(n0):
            raise FileNotFoundError(f"No N0.tsv found. Looked for:\n  {n0}")
    print(f"[hog] {n0}")

    # denominator HOGs and the target's status in each
    with open(n0, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        missing = [s for s in allsp if s not in header]
        if missing:
            sys.exit("Species not in N0.tsv: " + ", ".join(missing)
                     + "\nColumns: " + ", ".join(header[3:]))
        idx = {s: header.index(s) for s in allsp}

        n_total = 0
        denom_hogs = []          # rows passing the comparison-species m filter
        allpresent = set()       # target gene -> every comparison species present (any chrom)
        tgt_single_allpresent = set()
        unmapped = Counter()

        for row in reader:
            n_total += 1
            genes = {}
            for s in allsp:
                i = idx[s]
                cell = row[i].strip() if i < len(row) else ""
                genes[s] = [norm(x) for x in cell.split(",") if x.strip()]

            comp_present = all(genes[s] for s in comps)
            if comp_present and len(genes[target]) == 1:
                tgt_single_allpresent.add(genes[target][0])

            # comparison species: how many have >=1 gene on their m-chromosome
            n_on = 0
            comp_chroms = {}
            for s in comps:
                on_m, chroms_here = False, []
                for g in genes[s]:
                    c = pos[s].get(g)
                    if c is None:
                        unmapped[s] += 1
                        continue
                    chroms_here.append(c[0])
                    if c[0] == micro[s]:
                        on_m = True
                comp_chroms[s] = chroms_here
                if on_m:
                    n_on += 1
            if n_on < min_comp:
                continue

            tg = genes[target]
            tchroms = []
            for g in tg:
                c = pos[target].get(g)
                if c is None:
                    unmapped[target] += 1
                else:
                    tchroms.append(c[0])

            if not tg:
                status = "absent"
            elif len(tg) == 1:
                status = "micro_single" if tchroms == [tmicro] else "elsewhere_single"
            else:
                status = "micro_multi" if tmicro in tchroms else "elsewhere_multi"

            denom_hogs.append({"hog": row[0], "genes": genes, "status": status,
                               "tchroms": tchroms, "comp_chroms": comp_chroms,
                               "n_on": n_on})

    # counts and proportions
    for s in allsp:
        if unmapped[s]:
            print(f"[warn] {s}: {unmapped[s]} HOG genes had no BED entry (ID mismatch? --strip-re)")

    denom = len(denom_hogs)
    st = Counter(r["status"] for r in denom_hogs)
    numer = st["micro_single"]
    any_micro = st["micro_single"] + st["micro_multi"]

    print(f"\n[hog] total HOGs {n_total}")
    print(f"[hog] denominator (>= {min_comp} of {len(comps)} comparison species have >=1 gene "
          f"on their m): {denom}")
    if min_comp < len(comps):
        by_k = Counter(r["n_on"] for r in denom_hogs)
        print("[hog]   by number of comparison species on their m: "
              + ", ".join(f"{k}: {by_k[k]}" for k in sorted(by_k, reverse=True)))
    print(f"[hog] target status within denominator:")
    for k in ("micro_single", "micro_multi", "elsewhere_single", "elsewhere_multi", "absent"):
        print(f"      {k:17s} {st.get(k, 0):6d}")
    if denom:
        print(f"\n[proportion] single-copy on {tmicro}: {numer}/{denom} = {numer/denom:.4f}")
        print(f"[proportion] any copy on {tmicro}:     {any_micro}/{denom} = {any_micro/denom:.4f}")
        print(f"[proportion] present anywhere:        {denom - st.get('absent',0)}/{denom} "
              f"= {(denom - st.get('absent',0))/denom:.4f}")
    else:
        sys.exit("Denominator is 0: no HOG has every comparison species on its m-chromosome. "
                 "Check the m-chromosome IDs against each species' BED.")

    # where the target genes that are not on its m sit
    elsewhere = Counter(c for r in denom_hogs if r["status"].startswith("elsewhere")
                        for c in r["tchroms"])
    if elsewhere:
        print(f"\n[where] target chromosomes carrying denominator-HOG genes NOT on {tmicro}:")
        for c, n in elsewhere.most_common(8):
            print(f"      {c}: {n}")

    # shared flag: micro genes in the numerator; other classes single-copy with every comparison species present
    numer_genes = {r["genes"][target][0] for r in denom_hogs if r["status"] == "micro_single"}
    shared = {}
    for g, k in gclass.items():
        if k == "micro":
            if g in numer_genes:
                shared[g] = 1
        elif g in tgt_single_allpresent:
            shared[g] = 1

    # gene class table and denominator HOG table
    with open(args.out, "w") as fh:
        fh.write("gene\tchrom\tstart\tclass\tshared\n")
        for g, k in sorted(gclass.items()):
            c, s = pos[target][g]
            fh.write(f"{g}\t{c}\t{s}\t{k}\t{shared.get(g, 0)}\n")

    with open(args.hog_table, "w") as fh:
        cols = ["HOG", "status", "n_comp_on_m", f"{target}_genes", f"{target}_chroms"]
        for s in comps:
            cols += [f"{s}_genes", f"{s}_chroms"]
        fh.write("\t".join(cols) + "\n")
        for r in denom_hogs:
            out = [r["hog"], r["status"], str(r["n_on"]),
                   ";".join(r["genes"][target]), ";".join(r["tchroms"])]
            for s in comps:
                out += [";".join(r["genes"][s]), ";".join(r["comp_chroms"][s])]
            fh.write("\t".join(out) + "\n")

    # summary
    print("\n[sets] shared=1 genes per class (these feed the depth test):")
    cnt = Counter(gclass[g] for g in shared)
    tot = Counter(gclass.values())
    for k in CLS_ORDER:
        if k in tot:
            print(f"    {k:9s} {cnt.get(k, 0):6d} shared / {tot[k]:6d} total")
    if cnt.get("micro", 0) < 20:
        print(f"\n[warn] only {cnt.get('micro', 0)} micro genes in the numerator; "
              "report the CI, not a point value")
    print(f"\nwrote {args.out} and {args.hog_table}")


if __name__ == "__main__":
    main()
