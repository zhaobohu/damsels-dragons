#!/usr/bin/env python3
"""
CDS BED (chrom, start0, end, geneID) from a GFF3, each CDS attributed to its top-level ancestor (gene)
or, with --id-level transcript, to its direct parent. --check-bed reports ID overlap with a GENESPACE BED.

Example:
  python make_cds_bed.py --gff /path/to/Aeshna_juncea.annevo.gff3 \\
      --check-bed /path/to/palaeoptera/genespace/bed/Aeshna_juncea.bed \\
      --strip-re '\\.t\\d+$' --out cds.annevo.bed
"""
import argparse
import gzip
import re
import sys
from collections import Counter


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True)
    p.add_argument("--id-level", choices=["gene", "transcript"], default="gene",
                   help="attribute CDS to the top-level feature (gene) or its direct parent (transcript)")
    p.add_argument("--strip-re", default=None, help="regex stripped from the emitted IDs")
    p.add_argument("--check-bed", default=None,
                   help="GENESPACE BED for the same species; reports ID overlap")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    strip_re = re.compile(args.strip_re) if args.strip_re else None
    norm = (lambda g: strip_re.sub("", g)) if strip_re else (lambda g: g)

    # read the GFF: ID -> first Parent for every feature, and every CDS with its first Parent
    opener = gzip.open if args.gff.endswith(".gz") else open
    parent, cds = {}, []
    with opener(args.gff, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = {}
            for kv in f[8].strip().split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    a[k.strip()] = v.strip()
            if "ID" in a and "Parent" in a:
                parent[a["ID"]] = a["Parent"].split(",")[0]
            if f[2] == "CDS":
                pid = a.get("Parent")
                if pid is None:
                    continue
                cds.append((f[0], int(f[3]) - 1, int(f[4]), pid.split(",")[0]))

    if not cds:
        sys.exit("no CDS features with a Parent= attribute found -- is this GFF3? "
                 "(for GTF, convert first: gffread in.gtf -o out.gff3)")

    # attribute each CDS to its gene by walking Parent links up (guarding against loops), then sort
    rows, ids = [], set()
    for c, s, e, pid in cds:
        if args.id_level != "transcript":
            seen = set()
            while pid in parent and pid not in seen:
                seen.add(pid)
                pid = parent[pid]
        g = norm(pid)
        ids.add(g)
        rows.append((c, s, e, g))
    rows.sort(key=lambda r: (r[0], r[1]))

    # write the BED and a short summary
    with open(args.out, "w") as fh:
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")

    chroms = Counter(r[0] for r in rows)
    print(f"[cds] {len(rows)} CDS intervals, {len(ids)} genes -> {args.out}")
    print("[cds] sample IDs: " + ", ".join(sorted(ids)[:4]))
    print("[cds] top chromosomes: " + ", ".join(f"{c}({n})" for c, n in chroms.most_common(6)))

    # optional: ID overlap with the GENESPACE BED
    if args.check_bed:
        bed_ids = set()
        with open(args.check_bed) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    bed_ids.add(norm(f[3].strip()))
        hit = len(ids & bed_ids)
        print(f"[check] IDs shared with {args.check_bed}: {hit}/{len(ids)} "
              f"(BED has {len(bed_ids)})")
        if hit == 0:
            print("[check] ZERO overlap. Sample BED IDs: " + ", ".join(sorted(bed_ids)[:4]))
            print("        Try --id-level transcript, or a different --strip-re.")
            sys.exit(1)
        elif hit < 0.9 * len(ids):
            print("[check] WARNING: fewer than 90% of GFF genes match the BED; "
                  "check which side is dropping genes")


if __name__ == "__main__":
    main()
