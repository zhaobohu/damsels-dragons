#!/usr/bin/env python3
"""
Per-hit table from a RepeatMasker .align file with "Kimura (with divCpGMod)" lines (.gz read transparently).

    parse_rm_align.py genome.fa.align --chrom OZ369553.1 > repeats_chr12.tsv

Output columns: chrom start end strand repeat class family kimura id score
(start/end 1-based inclusive; kimura NA for blocks without a Kimura line).
"""
import argparse, gzip, re, sys

# alignment header line, and the Kimura line that follows it
HEAD = re.compile(r'^\s*\d+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+\S+\s+\d+\s+\d+\s+\(\d+\)')
KIM  = re.compile(r'^Kimura \(with divCpGMod\)\s*=\s*([\d.]+)')

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('align')
    ap.add_argument('--chrom', action='append', help='keep only these sequences (repeatable)')
    ap.add_argument('--min-len', type=int, default=0, help='drop hits shorter than this (bp)')
    a = ap.parse_args()
    keep = set(a.chrom) if a.chrom else None

    out = sys.stdout
    out.write('chrom\tstart\tend\tstrand\trepeat\tclass\tfamily\tkimura\tid\tscore\n')
    cur = None                 # header awaiting its Kimura line
    n_hit = n_kim = n_kept = 0

    # write the pending hit with its Kimura value, if it passes the filters
    def flush(kim):
        nonlocal cur, n_kept
        if cur is None:
            return
        chrom, s, e, strand, rep, cls, fam, rid, score = cur
        if (keep is None or chrom in keep) and (e - s + 1) >= a.min_len:
            out.write(f'{chrom}\t{s}\t{e}\t{strand}\t{rep}\t{cls}\t{fam}\t{kim}\t{rid}\t{score}\n')
            n_kept += 1
        cur = None

    # stream the .align, keeping only header lines and the Kimura line after each
    opener = gzip.open if a.align.endswith('.gz') else open
    with opener(a.align, 'rt') as fh:
        for line in fh:
            if HEAD.match(line):
                flush('NA')            # previous block had no Kimura line
                # score %div %del %ins query qbeg qend (qleft) [C] repeat#class/family r1 r2 (r3) id [n]
                # (the C strand flag shifts the repeat fields by one)
                t = line.split()
                score = t[0]
                chrom, qbeg, qend = t[4], int(t[5]), int(t[6])
                i = 8
                strand = '+'
                if t[i] == 'C':
                    strand = '-'; i += 1
                name = t[i]
                rid = t[i + 4] if len(t) > i + 4 else 'NA'
                if '#' in name:
                    rep, cf = name.split('#', 1)
                else:
                    rep, cf = name, 'Unknown'
                if '/' in cf:
                    cls, fam = cf.split('/', 1)
                else:
                    cls, fam = cf, ''
                cur = (chrom, qbeg, qend, strand, rep, cls, fam, rid, score)
                n_hit += 1
            elif cur is not None:
                m = KIM.match(line)
                if m:
                    n_kim += 1
                    flush(m.group(1))
        flush('NA')

    # summary to stderr
    sys.stderr.write(f'alignment blocks: {n_hit}   with Kimura: {n_kim}   written: {n_kept}\n')
    if n_hit and n_kim == 0:
        sys.stderr.write('WARNING: no Kimura lines found. Run calcDivergenceFromAlign.pl with -a first.\n')

if __name__ == '__main__':
    main()
