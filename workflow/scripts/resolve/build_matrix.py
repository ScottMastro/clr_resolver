#!/usr/bin/env python3
"""Collapse per-read allele calls into a per-barcode marker matrix.

A barcode marks one physical molecule = one orange copy. So we collapse all
reads of a barcode to a single marker vector:

  - SNV markers: the molecule's allele at a bubble = majority vote of its
    reads' calls (the ~5% per-read miscall averages out). Ties / unobserved
    sites are missing ('.').
  - Structural markers: how many of the molecule's read-path nodes fall in
    each orange sub-block a..e. A copy that lacks sub-block c never produces
    reads on orange_c nodes.

Outputs:
  --out-matrix  barcode, n_reads, n_called, then one column per informative
                bubble (allele 0/1, '.' missing)
  --out-struct  barcode, n_reads, covA..covE  (read-node counts per sub-block)
Both keep a `source_hap` column for convenience (NOT used by clustering).
"""

import argparse
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reads", required=True)
    ap.add_argument("--read-markers", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--node-color", required=True)
    ap.add_argument("--color", default="O")
    ap.add_argument("--out-matrix", required=True)
    ap.add_argument("--out-struct", required=True)
    ap.add_argument("--min-agree", type=float, default=0.7,
                    help="majority fraction needed to call a barcode allele")
    args = ap.parse_args()

    # informative bubbles, in order
    inf = []
    with open(args.markers) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[10] == "1":
                inf.append(c[0])
    inf_set = set(inf)
    print(f"[info] {len(inf)} informative markers", file=sys.stderr)

    # node -> orange sub-block letter
    sub = {}
    with open(args.node_color) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[1] == args.color and c[2].startswith("orange_"):
                sub[c[0]] = c[2].split("_")[1]

    # barcode -> source hap, read count, sub-block node counts
    bc_src = {}
    bc_reads = defaultdict(set)
    bc_cov = defaultdict(lambda: defaultdict(int))
    with open(args.reads) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            bc, src, rid, mate = c[8], c[2], c[0], c[1]
            bc_src[bc] = src
            bc_reads[bc].add((rid, mate))
            for n in c[11].split(","):
                if n in sub:
                    bc_cov[bc][sub[n]] += 1

    # barcode x bubble -> list of allele calls
    bc_bub = defaultdict(lambda: defaultdict(list))
    with open(args.read_markers) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rid, mate, bc, src, bid, al = c
            if bid in inf_set:
                bc_bub[bc][bid].append(int(al))

    barcodes = sorted(bc_src)

    # matrix
    with open(args.out_matrix, "w") as out:
        out.write("barcode\tsource_hap\tn_reads\tn_called\t"
                  + "\t".join(f"b{b}" for b in inf) + "\n")
        for bc in barcodes:
            row = []
            n_called = 0
            for b in inf:
                calls = bc_bub[bc].get(b, [])
                if not calls:
                    row.append(".")
                    continue
                ones = sum(calls)
                frac1 = ones / len(calls)
                if frac1 >= args.min_agree:
                    row.append("1")
                    n_called += 1
                elif frac1 <= 1 - args.min_agree:
                    row.append("0")
                    n_called += 1
                else:
                    row.append(".")
            out.write(f"{bc}\t{bc_src[bc]}\t{len(bc_reads[bc])}\t{n_called}\t"
                      + "\t".join(row) + "\n")

    # structural
    with open(args.out_struct, "w") as out:
        out.write("barcode\tsource_hap\tn_reads\tcovA\tcovB\tcovC\tcovD\tcovE\n")
        for bc in barcodes:
            cov = bc_cov[bc]
            out.write(f"{bc}\t{bc_src[bc]}\t{len(bc_reads[bc])}\t"
                      + "\t".join(str(cov.get(x, 0)) for x in "abcde") + "\n")

    print(f"[OK] {args.out_matrix}: {len(barcodes)} barcodes x {len(inf)} "
          f"markers", file=sys.stderr)
    print(f"[OK] {args.out_struct}", file=sys.stderr)


if __name__ == "__main__":
    main()
