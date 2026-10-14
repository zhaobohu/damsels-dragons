#!/usr/bin/env python3
"""GUIDANCE2-like column confidence scores for a pairwise PRANK codon alignment."""

import os
import sys
import subprocess
import tempfile
import shutil
import random
import argparse


def parse_fasta(path):
    seqs = []
    name, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(parts)))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        seqs.append((name, "".join(parts)))
    return seqs


def write_fasta(seqs, path):
    with open(path, "w") as f:
        for name, seq in seqs:
            f.write(f">{name}\n{seq}\n")


def run_prank(seq_file, out_prefix):
    cmd = ["prank", f"-d={seq_file}", f"-o={out_prefix}", "-codon", "-F", "-once"]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    out = f"{out_prefix}.best.fas"
    if not os.path.exists(out):
        raise FileNotFoundError(f"PRANK output not found: {out}")
    return out


def to_codons(seq):
    """Split an aligned sequence into codon triplets (gaps become '---')."""
    return [seq[i : i + 3] for i in range(0, len(seq), 3)]


def strip_gap_codons(codons):
    """Drop gap codons and return the ungapped nucleotide string."""
    return "".join(c for c in codons if c.replace("-", "") != "")


def residue_pairs(codons1, codons2):
    """(codon index in seq1, codon index in seq2) for every column where both are non-gap."""
    pairs = []
    ri, rj = 0, 0
    for c1, c2 in zip(codons1, codons2):
        g1 = c1.replace("-", "") == ""
        g2 = c2.replace("-", "") == ""
        if not g1 and not g2:
            pairs.append((ri, rj))
        if not g1:
            ri += 1
        if not g2:
            rj += 1
    return pairs


def main():
    ap = argparse.ArgumentParser(
        description="GUIDANCE2-like pairwise codon alignment confidence"
    )
    ap.add_argument("--seqfile", required=True, help="Input FASTA (2 unaligned CDS)")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--bootstraps", type=int, default=100)
    ap.add_argument("--col_cutoff", type=float, default=0.93)
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    os.makedirs(args.outdir, exist_ok=True)
    tmpdir = tempfile.mkdtemp(dir=args.outdir, prefix="bs_tmp_")

    try:
        # base PRANK codon alignment and its aligned codon pairs
        base_prefix = os.path.join(args.outdir, "base")
        base_aln = run_prank(args.seqfile, base_prefix)
        seqs = parse_fasta(base_aln)
        names = [s[0] for s in seqs]
        base_c1 = to_codons(seqs[0][1])
        base_c2 = to_codons(seqs[1][1])
        n_codons = len(base_c1)

        base_pairs = residue_pairs(base_c1, base_c2)
        n_paired = len(base_pairs)

        if n_paired == 0:
            sys.stderr.write("No aligned codon pairs in base alignment.\n")
            print("NA\tNA")
            return

        # bootstrap: resample codon columns, strip gaps, realign, count reproduced pairs
        hits = [0] * n_paired
        completed = 0

        for b in range(args.bootstraps):
            idx = random.choices(range(n_codons), k=n_codons)
            bs_c1 = [base_c1[i] for i in idx]
            bs_c2 = [base_c2[i] for i in idx]

            bs_seq1 = strip_gap_codons(bs_c1)
            bs_seq2 = strip_gap_codons(bs_c2)

            # need at least 2 codons per sequence
            if len(bs_seq1) < 6 or len(bs_seq2) < 6:
                continue

            bs_fa = os.path.join(tmpdir, f"bs_{b}.fa")
            bs_prefix = os.path.join(tmpdir, f"bs_{b}")
            write_fasta([(names[0], bs_seq1), (names[1], bs_seq2)], bs_fa)

            # replicates where PRANK fails are skipped
            try:
                bs_aln = run_prank(bs_fa, bs_prefix)
                bs_seqs = parse_fasta(bs_aln)
                bs_pairs = set(residue_pairs(to_codons(bs_seqs[0][1]),
                                             to_codons(bs_seqs[1][1])))
            except Exception:
                continue

            for i, pair in enumerate(base_pairs):
                if pair in bs_pairs:
                    hits[i] += 1
            completed += 1

            # remove replicate files to save inodes
            for ext in [".fa", ".best.fas", ".best.dnd", ".best.xml",
                        ".1.fas", ".1.dnd", ".1.xml"]:
                p = bs_prefix + ext
                if os.path.exists(p):
                    os.remove(p)
            if os.path.exists(bs_fa):
                os.remove(bs_fa)

        if completed == 0:
            sys.stderr.write("All bootstrap replicates failed.\n")
            print("NA\tNA")
            return

        col_scores = [h / completed for h in hits]

        # per-column scores for the aligned codon pairs
        score_file = os.path.join(args.outdir, "column_scores.tsv")
        with open(score_file, "w") as f:
            f.write("aln_codon_col\tres_i\tres_j\tscore\n")
            pi = 0
            ri, rj = 0, 0
            for ci, (c1, c2) in enumerate(zip(base_c1, base_c2)):
                g1 = c1.replace("-", "") == ""
                g2 = c2.replace("-", "") == ""
                if not g1 and not g2:
                    f.write(f"{ci}\t{ri}\t{rj}\t{col_scores[pi]:.4f}\n")
                    pi += 1
                if not g1:
                    ri += 1
                if not g2:
                    rj += 1

        # filtered alignment: drop paired columns scoring below the cutoff, keep gap columns
        pi = 0
        keep_cols = []
        for ci in range(n_codons):
            c1 = base_c1[ci]
            c2 = base_c2[ci]
            g1 = c1.replace("-", "") == ""
            g2 = c2.replace("-", "") == ""
            if g1 or g2:
                keep_cols.append(ci)       # yn00 ignores gap columns anyway
            else:
                if col_scores[pi] >= args.col_cutoff:
                    keep_cols.append(ci)
                pi += 1

        filt_seq1 = "".join(base_c1[i] for i in keep_cols)
        filt_seq2 = "".join(base_c2[i] for i in keep_cols)
        filt_path = os.path.join(args.outdir, "filtered.fa")
        write_fasta([(names[0], filt_seq1), (names[1], filt_seq2)], filt_path)

        # summary table
        n_retained = sum(1 for s in col_scores if s >= args.col_cutoff)
        mean_score = sum(col_scores) / n_paired
        frac_retained = n_retained / n_paired

        with open(os.path.join(args.outdir, "summary.tsv"), "w") as f:
            f.write(f"mean_score\t{mean_score:.4f}\n")
            f.write(f"frac_retained\t{frac_retained:.4f}\n")
            f.write(f"n_paired_cols\t{n_paired}\n")
            f.write(f"n_retained_cols\t{n_retained}\n")
            f.write(f"bootstraps_completed\t{completed}\n")
            f.write(f"col_cutoff\t{args.col_cutoff}\n")

        # mean score and retained fraction on stdout for the calling shell script
        print(f"{mean_score:.4f}\t{frac_retained:.4f}")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
