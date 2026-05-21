#!/usr/bin/env python3
"""Per-CLR-read truth: which orange copy did the read come from.

Two inputs give the answer:
  - the pbsim `.maf` per haplotype: each block's `s ref <start> <len> ...`
    line is the read's true span in that haplotype's window coordinates.
  - the haplotype's `merged.tsv` block decomposition: orange_* rows in the
    same window coordinates. Coordinate-adjacent orange_* rows form one
    orange copy; the copy's signature is the set of sub-blocks it carries.

A read is assigned to the copy its true span overlaps most. Copies are
numbered c1..cN per haplotype in coordinate order.

Output (TSV, score.py-compatible — col0 read id, col4 truth copy):
  read_id, source_hap, n_aln, mol_span, truth_copy, copy_sig
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_maf(path):
    """read name (S1_<n>) -> (ref_start, ref_end) of its longest block."""
    spans = {}
    with open_maybe_gz(path) as fh:
        ref = None
        for line in fh:
            if line.startswith("s "):
                f = line.split()
                # f = ['s', name, start, len, strand, srcsize, seq]
                name, start, length = f[1], int(f[2]), int(f[3])
                if name == "ref":
                    ref = (start, start + length)
                elif ref is not None:
                    s, e = ref
                    if name not in spans or (e - s) > (spans[name][1] -
                                                       spans[name][0]):
                        spans[name] = (s, e)
                    ref = None
    return spans


def copies_from_paint(gfa, node_color_path, sample, color, gap, min_copy):
    """Orange-copy intervals straight from the graph painting.

    The matrix is built from the node-colour painting (node_color.tsv, itself
    a projection of hg38.haplogroup.bed onto the graph). merged.tsv is only an
    approximate per-haplotype redrawing whose block boundaries drift ~8 kb from
    the painting -- so truth must use the SAME painting the matrix does.

    Each haplotype is a P-line walk through the graph; accumulating node
    lengths gives window coordinates, and maximal runs of `color`-painted
    nodes (small gaps bridged) are the orange copies.

    Returns {hap_tag: [(start, end), ...]} in window coordinates.
    """
    col = {}
    with open(node_color_path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            col[c[0]] = c[1]
    nlen, plines = {}, {}
    with open(gfa) as fh:
        for line in fh:
            if line[0] == "S":
                p = line.rstrip("\n").split("\t")
                nlen[p[1]] = len(p[2]) if len(p) > 2 and p[2] != "*" else 0
            elif line[0] == "P":
                p = line.split("\t")
                if p[1].startswith(sample):
                    plines[p[1]] = p[2]
    out = {}
    for name, path in plines.items():
        tag = "h1" if "#1#" in name else ("h2" if "#2#" in name else None)
        if tag is None:
            continue
        off, raw = 0, []
        for step in path.split(","):
            n = step[:-1]
            L = nlen.get(n, 0)
            if col.get(n) == color and L > 0:
                raw.append((off, off + L))
            off += L
        merged = []
        for s, e in raw:
            if merged and s - merged[-1][1] <= gap:
                merged[-1][1] = e
            else:
                merged.append([s, e])
        out[tag] = [(s, e) for s, e in merged if e - s >= min_copy]
    return out


def load_copies(merged, copy_gap):
    """merged.tsv -> list of (start, end, signature) orange copies in order."""
    oranges = []
    with open(merged) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 4 or not c[1].startswith("orange"):
                continue
            oranges.append((int(c[2]), int(c[3]), c[1].split("_")[-1]))
    oranges.sort()
    if not oranges:
        return []
    runs, cur = [], [oranges[0]]
    for row in oranges[1:]:
        if row[0] - cur[-1][1] <= copy_gap:
            cur.append(row)
        else:
            runs.append(cur)
            cur = [row]
    runs.append(cur)
    return [(run[0][0], run[-1][1],
             "".join(sorted(set(sub for _, _, sub in run)))) for run in runs]


def overlap(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maf-h1", required=True)
    ap.add_argument("--maf-h2", required=True)
    ap.add_argument("--gfa", required=True,
                    help="graph GFA with per-haplotype P-lines")
    ap.add_argument("--node-color", required=True,
                    help="node_color.tsv -- the painting the matrix uses")
    ap.add_argument("--sample", default="HG00133",
                    help="P-line name prefix of the sample")
    ap.add_argument("--color", default="O")
    ap.add_argument("--out", required=True)
    ap.add_argument("--paint-gap", type=int, default=2000,
                    help="bridge orange runs separated by < this many bp")
    ap.add_argument("--min-copy", type=int, default=3000,
                    help="drop painted orange runs shorter than this")
    ap.add_argument("--merged-h1",
                    help="optional: merged.tsv, used only for the copy_sig "
                         "diagnostic column (order/content, not boundaries)")
    ap.add_argument("--merged-h2")
    ap.add_argument("--copy-gap", type=int, default=20000)
    args = ap.parse_args()

    # copy boundaries: from the graph painting (NOT merged.tsv -- its block
    # boundaries drift ~8 kb from the painting the matrix is built on).
    paint = copies_from_paint(args.gfa, args.node_color, args.sample,
                              args.color, args.paint_gap, args.min_copy)
    mafs = {"h1": load_maf(args.maf_h1), "h2": load_maf(args.maf_h2)}

    # copy_sig is diagnostic only -- pull it from merged.tsv by copy order.
    sigs = {"h1": [], "h2": []}
    for tag, mp in (("h1", args.merged_h1), ("h2", args.merged_h2)):
        if mp:
            sigs[tag] = [s for _, _, s in load_copies(mp, args.copy_gap)]

    for tag in ("h1", "h2"):
        cps = paint.get(tag, [])
        desc = " ".join(f"c{i+1}:{c0}-{c1}" for i, (c0, c1) in enumerate(cps))
        print(f"[info] {tag}: {len(mafs[tag])} reads, painted copies = {desc}",
              file=sys.stderr)

    n_assigned = defaultdict(int)
    with open(args.out, "w") as out:
        out.write("read_id\tsource_hap\tn_aln\tmol_span\ttruth_copy\t"
                  "copy_sig\n")
        for tag in ("h1", "h2"):
            copies = paint.get(tag, [])
            for name, (s, e) in sorted(mafs[tag].items()):
                rid = f"{tag}_{name}"
                best_i, best_ov = -1, 0
                for i, (c0, c1) in enumerate(copies):
                    ov = overlap(s, e, c0, c1)
                    if ov > best_ov:
                        best_i, best_ov = i, ov
                if best_i < 0:
                    tc, csig = ".", "."
                else:
                    tc = f"{tag}#c{best_i+1}"
                    csig = (sigs[tag][best_i]
                            if best_i < len(sigs[tag]) else ".")
                    n_assigned[tc] += 1
                out.write(f"{rid}\t{tag}\t1\t{s}:{e}\t{tc}\t{csig}\n")

    print(f"[OK] {args.out}: copy sizes = {dict(sorted(n_assigned.items()))}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
