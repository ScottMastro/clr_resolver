#!/usr/bin/env python3
"""Infer the orange (MUC20) copy number of a CLR read set by cosigt-style
coverage matching against the graph's panel haplotypes.

A sample is diploid: its reads are explained by some pair of panel haplotypes
(the pair may be one haplotype taken twice -- a homozygous match). Each panel
haplotype's orange segment is a known copy signature: the set of orange nodes
its P-line walks. This script takes the sample's per-orange-node read depth
(from a GraphAligner GAF -- no separate coverage tool needed), finds the
haplotype pair whose summed signature best explains it (cosine similarity),
and reports the copy number k = orange copies(hap1) + orange copies(hap2),
where a haplotype's orange copy count is its painted-run count.

Anchor nodes set the depth scale: their median depth is two haplotype
traversals, so per-haplotype depth d = anchor_depth / 2; the matched pair's
predicted orange coverage is cross-checked against the observed total at d.

For an honest novel-sample estimate, --exclude-sample drops the sample's own
panel haplotypes (leave-one-out).

Output schema (--out): metric, value -- rows estimated_k, matched_hap1,
matched_hap2, copies_hap1, copies_hap2, match_cosine, anchor_depth,
depth_per_haplotype, observed_over_predicted."""

import argparse
import os
import re
import sys

import numpy as np

_SEG = re.compile(r"(\d+)[+-]")
_STEP = re.compile(r"[<>](\d+)")


def load_painting(path):
    """node id -> base colour."""
    col = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2:
                col[c[0]] = c[1]
    return col


def parse_graph(gfa):
    """Return (node_len, paths) where paths is {path_name: [node ids]}."""
    node_len, paths = {}, {}
    with open(gfa) as fh:
        for line in fh:
            t = line[0]
            if t == "S":
                f = line.split("\t")
                node_len[f[1]] = len(f[2].rstrip("\n")) if len(f) > 2 else 0
            elif t == "P":
                f = line.split("\t")
                paths[f[1]] = _SEG.findall(f[2])
            elif t == "W":
                f = line.split("\t")
                name = "#".join(f[1:4])
                paths[name] = _STEP.findall(f[6])
    return node_len, paths


def orange_runs(path, node_len, col, color, gap=2000, min_copy=3000):
    """Count painted copies of `color` along a path (bridge < gap, drop short)."""
    off, runs = 0, []
    for n in path:
        L = node_len.get(n, 0)
        if col.get(n, ".") == color:
            if runs and off - runs[-1][1] < gap:
                runs[-1][1] = off + L
            else:
                runs.append([off, off + L])
        off += L
    return sum(1 for s, e in runs if e - s >= min_copy)


def gaf_depth(gaf):
    """node id -> number of GAF alignments crossing it."""
    depth = {}
    with open(gaf) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 6:
                continue
            for n in _STEP.findall(f[5]):
                depth[n] = depth.get(n, 0) + 1
    return depth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaf", required=True, help="GraphAligner GAF of the reads")
    ap.add_argument("--gfa", required=True, help="graph GFA (panel haplotypes)")
    ap.add_argument("--node-color", required=True, help="node_color.tsv")
    ap.add_argument("--anchor-nodes", required=True,
                    help="anchor-node TSV from find_anchor_nodes.py")
    ap.add_argument("--color", default="O", help="target colour block")
    ap.add_argument("--exclude-sample", default=None,
                    help="drop panel haplotypes of this sample (leave-one-out)")
    ap.add_argument("--out", required=True, help="metric/value result TSV")
    args = ap.parse_args()

    col = load_painting(args.node_color)
    node_len, paths = parse_graph(args.gfa)
    depth = gaf_depth(args.gaf)

    anchors = []
    with open(args.anchor_nodes) as fh:
        fh.readline()
        for line in fh:
            anchors.append(line.split("\t")[0])
    anchor_d = np.median([depth.get(n, 0) for n in anchors]) if anchors else 0.0
    d = anchor_d / 2.0                       # per-haplotype depth

    color_nodes = sorted((n for n in node_len if col.get(n, ".") == args.color),
                         key=int)
    idx = {n: i for i, n in enumerate(color_nodes)}
    if not color_nodes:
        sys.exit(f"[ERROR] estimate_copy_number: no '{args.color}' nodes")

    # panel: one orange signature + copy count per haplotype (sample dropped)
    names, sig, copies = [], [], []
    for name, path in paths.items():
        if args.exclude_sample and name.startswith(args.exclude_sample + "#"):
            continue
        v = np.zeros(len(color_nodes), dtype=np.float64)
        for n in path:
            if n in idx:
                v[idx[n]] += 1.0
        if v.sum() == 0:
            continue
        names.append(name)
        sig.append(v)
        copies.append(orange_runs(path, node_len, col, args.color))
    if not names:
        sys.exit("[ERROR] estimate_copy_number: no panel haplotypes left")
    SIG = np.vstack(sig)                                  # H x nodes
    copies = np.array(copies)

    obs = np.array([depth.get(n, 0) for n in color_nodes], dtype=np.float64)
    if d <= 0:
        sys.exit("[ERROR] estimate_copy_number: anchor depth is zero")

    # Match in depth-calibrated units: obs/d is the per-node copy count, and a
    # haplotype pair predicts sig_a + sig_b. Pick the pair minimising
    # ||obs/d - (sig_a + sig_b)||^2 -- this keeps the depth magnitude (cosine
    # would discard it and the copy number with it).
    t = obs / d
    th = SIG @ t                                          # H
    G = SIG @ SIG.T                                       # H x H
    diag = np.diag(G)
    # ||t - sig_a - sig_b||^2  =  |t|^2 - 2(th_a+th_b) + diag_a+diag_b+2 G_ab
    err = (-2 * (th[:, None] + th[None, :])
           + diag[:, None] + diag[None, :] + 2 * G)
    a, b = np.unravel_index(np.argmin(err), err.shape)

    k = int(copies[a] + copies[b])
    pred_sum = d * float((SIG[a] + SIG[b]).sum())
    obs_pred = obs.sum() / pred_sum if pred_sum > 0 else 0.0

    # depth-only cross-check: total orange traversals / d / nodes-per-copy
    per_copy = np.median([SIG[i].sum() / copies[i]
                          for i in range(len(names)) if copies[i] > 0])
    k_depth = obs.sum() / d / per_copy if per_copy > 0 else 0.0

    rows = [
        ("estimated_k", k),
        ("matched_hap1", names[a]),
        ("matched_hap2", names[b]),
        ("copies_hap1", int(copies[a])),
        ("copies_hap2", int(copies[b])),
        ("match_resid_rms", f"{np.sqrt(max(err[a, b], 0) / len(t)):.4f}"),
        ("anchor_depth", f"{anchor_d:.1f}"),
        ("depth_per_haplotype", f"{d:.1f}"),
        ("observed_over_predicted", f"{obs_pred:.3f}"),
        ("k_depth_crosscheck", f"{k_depth:.2f}"),
    ]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as out:
        out.write("metric\tvalue\n")
        for m, v in rows:
            out.write(f"{m}\t{v}\n")

    sys.stderr.write(
        f"[INFO] estimate_copy_number: matched {names[a]} (k={copies[a]}) + "
        f"{names[b]} (k={copies[b]}), obs/pred={obs_pred:.2f}, "
        f"k_depth={k_depth:.2f}\n")
    sys.stderr.write(f"[OK] estimate_copy_number: estimated_k={k} -> {args.out}\n")


if __name__ == "__main__":
    main()
