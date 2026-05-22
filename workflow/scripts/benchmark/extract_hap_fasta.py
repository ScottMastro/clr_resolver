#!/usr/bin/env python3
"""Extract a sample's two haplotype region assemblies from a pangenome graph.

Each haplotype is a P-line walk through the graph; concatenating the segment
sequences in path order (reverse-complemented on '-' steps) reconstructs that
haplotype's region sequence. Used to give the benchmark a known assembly to
simulate reads from when no standalone FASTA is on hand.

P-lines are matched by the `<sample>#<1|2>#` PanSN prefix.

Output: one FASTA per haplotype (--out-hap1, --out-hap2)."""

import argparse
import os
import sys

_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_segments(gfa):
    """node id -> sequence (S-lines)."""
    seq = {}
    with open(gfa) as fh:
        for line in fh:
            if line[0] == "S":
                f = line.split("\t")
                seq[f[1]] = f[2].rstrip("\n")
    return seq


def find_pline(gfa, prefix):
    """The first P-line whose name starts with prefix -> (name, steps)."""
    with open(gfa) as fh:
        for line in fh:
            if line[0] != "P":
                continue
            f = line.split("\t")
            if f[1].startswith(prefix):
                return f[1], f[2]
    return None, None


def build(seq, steps):
    """Concatenate segment sequences along a P-line step string."""
    out = []
    for step in steps.split(","):
        nid, orient = step[:-1], step[-1]
        s = seq.get(nid, "")
        if orient == "-":
            s = s.translate(_RC)[::-1]
        out.append(s)
    return "".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True, help="input graph GFA")
    ap.add_argument("--sample", required=True, help="sample id (PanSN prefix)")
    ap.add_argument("--out-hap1", required=True, help="haplotype-1 FASTA")
    ap.add_argument("--out-hap2", required=True, help="haplotype-2 FASTA")
    args = ap.parse_args()

    seq = parse_segments(args.gfa)
    for hap, out in ((1, args.out_hap1), (2, args.out_hap2)):
        name, steps = find_pline(args.gfa, f"{args.sample}#{hap}#")
        if name is None:
            sys.exit(f"[ERROR] extract_hap_fasta: no P-line "
                     f"{args.sample}#{hap}# in {args.gfa}")
        s = build(seq, steps)
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        with open(out, "w") as fh:
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
        sys.stderr.write(f"[INFO] extract_hap_fasta: {name} {len(s):,} bp\n")
    sys.stderr.write(f"[OK] extract_hap_fasta: {args.sample} -> "
                     f"{args.out_hap1}, {args.out_hap2}\n")


if __name__ == "__main__":
    main()
