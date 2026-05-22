#!/usr/bin/env python3
"""Find a pangenome graph's anchor nodes and summarise path traversal counts.

An anchor node is a segment present in (nearly) every path -- a well-covered
node whose coverage is defined for every haplotype, so it can serve as a basis
for coverage matching. Anchors are NOT restricted to single-copy: a node every
copy of a region traverses is visited once per copy, so its per-path
multiplicity carries the copy-number signal. Multiplicity is reported, not
filtered, so callers choose the role:

  - single-copy anchors (max_mult == 1) are the depth reference -- in a
    diploid sample their read depth is exactly two haplotype traversals;
  - multi-copy anchors (max_mult > 1) carry the copy-number signal.

Reads P-lines and W-lines from a GFA. A node is an anchor when it appears in
at least --min-path-frac of paths; --max-mult optionally caps multiplicity.

Output schema (--out): node, length, n_paths, path_frac, total_occ,
max_mult, mean_mult[, color]. A traversal-count summary is written to stderr."""

import argparse
import os
import re
import sys
from collections import Counter

_SEG = re.compile(r"(\d+)[+-]")          # P-line: 123+,124-,...
_STEP = re.compile(r"[<>](\d+)")         # W-line walk: >123<124...


def parse_graph(gfa):
    """Return (node_len, paths) where paths is a list of node-id lists."""
    node_len, paths = {}, []
    with open(gfa) as fh:
        for line in fh:
            t = line[0]
            if t == "S":
                f = line.split("\t")
                node_len[f[1]] = len(f[2].rstrip("\n")) if len(f) > 2 else 0
            elif t == "P":
                f = line.split("\t")
                paths.append(_SEG.findall(f[2]))
            elif t == "W":
                f = line.split("\t")
                paths.append(_STEP.findall(f[6]))
    return node_len, paths


def tally(paths):
    """Per node: (#paths containing it, total occurrences, max multiplicity)."""
    n_paths, total_occ, max_mult = Counter(), Counter(), {}
    for nodes in paths:
        seen = Counter(nodes)
        for n, m in seen.items():
            n_paths[n] += 1
            total_occ[n] += m
            if m > max_mult.get(n, 0):
                max_mult[n] = m
    return n_paths, total_occ, max_mult


def summary(node_len, npath, n_paths, max_mult, anchors, color):
    """Write a traversal-count summary to stderr."""
    w = sys.stderr.write
    w(f"[INFO] find_anchor_nodes: {len(node_len)} nodes, {npath} paths\n")

    buckets = [(1.0, "in every path"), (0.99, ">=99%"), (0.90, ">=90%"),
               (0.50, ">=50%"), (0.0, "<50%")]
    fracs = sorted((n_paths.get(n, 0) / npath for n in node_len), reverse=True)
    w("[INFO] node ubiquity:\n")
    lo = 1.01
    for thr, label in buckets:
        c = sum(1 for f in fracs if thr <= f < lo)
        w(f"[INFO]   {label:14s} {c:7d} nodes\n")
        lo = thr

    w(f"[INFO] anchor nodes: {len(anchors)}"
      f" ({sum(node_len.get(n, 0) for n in anchors):,} bp)\n")
    mm = Counter(max_mult.get(n, 0) for n in anchors)
    w("[INFO] anchor max-multiplicity: "
      + ", ".join(f"{k}x:{mm[k]}" for k in sorted(mm)) + "\n")
    if color:
        cc = Counter(color.get(n, ".") for n in anchors)
        w("[INFO] anchor colours: "
          + ", ".join(f"{k}:{cc[k]}" for k in sorted(cc)) + "\n")
        # colour x single/multi-copy cross-tab
        for c in sorted(cc):
            sub = [n for n in anchors if color.get(n, ".") == c]
            single = sum(1 for n in sub if max_mult.get(n, 0) == 1)
            multi = len(sub) - single
            hi = max((max_mult.get(n, 0) for n in sub), default=0)
            w(f"[INFO]   colour {c}: {single} single-copy, {multi} multi-copy"
              f" (max multiplicity {hi})\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True, help="input graph GFA")
    ap.add_argument("--out", required=True, help="anchor-node TSV")
    ap.add_argument("--min-path-frac", type=float, default=1.0,
                    help="min fraction of paths a node must appear in")
    ap.add_argument("--max-mult", type=int, default=0,
                    help="optional cap on per-path multiplicity (0 = no cap)")
    ap.add_argument("--node-color", default=None,
                    help="optional node_color.tsv to add a colour column")
    args = ap.parse_args()

    node_len, paths = parse_graph(args.gfa)
    if not paths:
        sys.exit("[ERROR] find_anchor_nodes: no P-lines or W-lines in GFA")
    npath = len(paths)
    n_paths, total_occ, max_mult = tally(paths)

    color = {}
    if args.node_color:
        with open(args.node_color) as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) >= 2:
                    color[c[0]] = c[1]

    anchors = sorted(
        (n for n in node_len
         if n_paths.get(n, 0) > 0
         and n_paths[n] / npath >= args.min_path_frac
         and (args.max_mult == 0 or max_mult.get(n, 0) <= args.max_mult)),
        key=int)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        hdr = "node\tlength\tn_paths\tpath_frac\ttotal_occ\tmax_mult\tmean_mult"
        out.write(hdr + ("\tcolor\n" if args.node_color else "\n"))
        for n in anchors:
            row = (f"{n}\t{node_len.get(n, 0)}\t{n_paths[n]}\t"
                   f"{n_paths[n] / npath:.4f}\t{total_occ[n]}\t"
                   f"{max_mult.get(n, 0)}\t{total_occ[n] / n_paths[n]:.3f}")
            if args.node_color:
                row += f"\t{color.get(n, '.')}"
            out.write(row + "\n")

    summary(node_len, npath, n_paths, max_mult, anchors, color)
    sys.stderr.write(f"[OK] find_anchor_nodes: wrote {args.out} "
                     f"({len(anchors)} anchor nodes)\n")


if __name__ == "__main__":
    main()
