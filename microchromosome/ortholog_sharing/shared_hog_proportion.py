#!/usr/bin/env python3
"""
Proportion of HOGs present in every --shared species that are also present in --target
(OrthoFinder N0.tsv; species not listed are ignored).

Example:
    python shared_hog_proportion.py \\
        --orthofinder-dir /path/to/OrthoFinder/Results_XXX \\
        --shared Pantala_flavescens Ischnura_elegans Platycnemis_pennipes \\
        --target Aeshna_juncea
"""
import argparse
import csv
import os
import sys


def cell_present(row, i):
    return i < len(row) and row[i].strip() != ""


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--orthofinder-dir", required=True,
                   help="OrthoFinder Results_* directory (or a direct path to N0.tsv)")
    p.add_argument("--shared", nargs="+", default=[],
                   help="Species columns defining the denominator set")
    p.add_argument("--target", default=None,
                   help="Extra species added for the numerator")
    p.add_argument("--list-species", action="store_true",
                   help="Print the species column names in N0.tsv and exit")
    args = p.parse_args()

    # N0.tsv: given directly, or inside the OrthoFinder Results_* directory
    if os.path.isfile(args.orthofinder_dir):
        n0 = args.orthofinder_dir
    else:
        n0 = os.path.join(args.orthofinder_dir, "Phylogenetic_Hierarchical_Orthogroups", "N0.tsv")
        if not os.path.isfile(n0):
            raise FileNotFoundError(
                f"Could not find N0.tsv. Looked for:\n  {n0}\n"
                "Pass --orthofinder-dir pointing at the OrthoFinder Results_* directory, "
                "or give the path to N0.tsv directly."
            )

    with open(n0, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        # species columns follow HOG, OG, Gene Tree Parent Clade
        species_cols = header[3:]

        # optionally just list the species columns
        if args.list_species:
            print("Species columns in N0.tsv:")
            for s in species_cols:
                print(f"  {s}")
            return

        # check the requested species are columns
        requested = list(args.shared) + ([args.target] if args.target else [])
        missing = [s for s in requested if s not in species_cols]
        if missing:
            sys.exit(
                "These requested species are not columns in N0.tsv:\n  "
                + "\n  ".join(missing)
                + "\n\nAvailable columns:\n  "
                + "\n  ".join(species_cols)
                + "\n\n(Column names come from your input FASTA filenames -- "
                  "match them exactly.)"
            )
        if not args.shared:
            sys.exit("Provide --shared (and usually --target). "
                     "Use --list-species to see the names.")

        # count HOGs present in all shared species, and of those, present in the target
        shared_idx = [header.index(s) for s in args.shared]
        target_idx = header.index(args.target) if args.target else None

        total = denom = numer = 0
        for row in reader:
            total += 1
            if all(cell_present(row, i) for i in shared_idx):
                denom += 1
                if target_idx is not None and cell_present(row, target_idx):
                    numer += 1

    # report
    print(f"N0.tsv: {n0}")
    print(f"Total HOGs: {total}")
    print(f"Shared species (denominator set): {', '.join(args.shared)}")
    if args.target:
        print(f"Target added for numerator: {args.target}")
    print()
    print(f"Denominator (HOGs shared by the {len(args.shared)} species): {denom}")
    if args.target:
        if denom == 0:
            print("Denominator is 0 -- cannot compute a ratio. Check species names.")
        else:
            ratio = numer / denom
            print(f"Numerator   (+ {args.target}): {numer}")
            print(f"Proportion  numerator / denominator: {ratio:.4f}  ({numer}/{denom})")


if __name__ == "__main__":
    main()
