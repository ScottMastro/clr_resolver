#!/usr/bin/env python3
"""Leave-one-out evaluation of orange copy-number estimation on the panel.

For every sample, its two haplotype P-lines are the noise-free "observation";
the sample is dropped from the panel and its orange copy number k is inferred
by cosine-matching the observation to the best remaining haplotype pair. k is
read off the matched pair's painted copy count. Reports per-k and overall
accuracy -- the graph-theoretic ceiling of the matcher before read noise.

The matching basis is the well-covered anchor set (nodes in every path),
multiplicity-valued; --max-mult drops hyper-repeat (VNTR) nodes.

Output: accuracy table to stdout."""

import argparse
import re
import sys
from collections import Counter

import numpy as np

_SEG = re.compile(r"(\d+)[+-]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--node-color", required=True)
    ap.add_argument("--max-mult", type=int, default=12,
                    help="drop anchor nodes whose max per-path multiplicity "
                         "exceeds this (VNTR filter)")
    ap.add_argument("--basis", choices=["all", "orange"], default="all",
                    help="match on all anchor nodes or orange anchors only")
    ap.add_argument("--metric", choices=["cosine", "lsq"], default="lsq",
                    help="cosine similarity (scale-free) or least-squares "
                         "(uses copy-count magnitude)")
    args = ap.parse_args()

    col, nlen = {}, {}
    for ln in open(args.node_color).readlines()[1:]:
        c = ln.rstrip("\n").split("\t")
        col[c[0]] = c[1]
    paths = {}
    for ln in open(args.gfa):
        if ln[0] == "S":
            p = ln.split("\t")
            nlen[p[1]] = len(p[2].rstrip("\n")) if len(p) > 2 else 0
        elif ln[0] == "P":
            p = ln.split("\t")
            paths[p[1]] = _SEG.findall(p[2])

    def ocopies(path, gap=2000, minc=3000):
        off, runs = 0, []
        for n in path:
            L = nlen.get(n, 0)
            if col.get(n, ".") == "O":
                if runs and off - runs[-1][1] < gap:
                    runs[-1][1] = off + L
                else:
                    runs.append([off, off + L])
            off += L
        return sum(1 for s, e in runs if e - s >= minc)

    names = [n for n in paths if not n.startswith(("GRCh38", "CHM13"))]
    cnt = {n: Counter(paths[n]) for n in names}
    npath = len(names)

    # anchor set: nodes in every path; VNTR filter; optional orange-only
    in_all = [n for n in nlen
              if sum(1 for nm in names if n in cnt[nm]) == npath]
    maxmult = {n: max(cnt[nm].get(n, 0) for nm in names) for n in in_all}
    basis = [n for n in in_all if maxmult[n] <= args.max_mult
             and (args.basis == "all" or col.get(n, ".") == "O")]
    idx = {n: i for i, n in enumerate(basis)}
    sys.stderr.write(f"[INFO] {npath} haplotypes, anchor basis {len(basis)} "
                     f"nodes (max-mult {args.max_mult}, {args.basis})\n")

    SIG = np.zeros((npath, len(basis)), dtype=np.float64)
    for i, nm in enumerate(names):
        for n, m in cnt[nm].items():
            if n in idx:
                SIG[i, idx[n]] = m
    ocop = np.array([ocopies(paths[nm]) for nm in names])

    # sample -> its two haplotype row indices
    samp = {}
    for i, nm in enumerate(names):
        samp.setdefault(nm.split("#")[0], []).append(i)
    samp = {s: ix for s, ix in samp.items() if len(ix) == 2}

    G = SIG @ SIG.T
    diag = np.diag(G)
    correct = 0
    by_k = Counter()
    by_k_ok = Counter()
    misses = []
    resids = []
    for s, (i1, i2) in sorted(samp.items()):
        obs = SIG[i1] + SIG[i2]
        k_true = int(ocop[i1] + ocop[i2])
        oh = SIG @ obs
        if args.metric == "cosine":
            # cosine(obs, sig_a+sig_b): scale-free, ignores copy magnitude
            num = oh[:, None] + oh[None, :]
            den = np.sqrt(np.maximum(
                diag[:, None] + diag[None, :] + 2 * G, 1e-9)) \
                * np.sqrt(obs @ obs)
            score = num / den
            score[i1, :] = score[:, i1] = -1
            score[i2, :] = score[:, i2] = -1
            a, b = np.unravel_index(np.argmax(score), score.shape)
        else:
            # ||obs - (sig_a+sig_b)||^2: keeps the copy-count magnitude
            err = (-2 * (oh[:, None] + oh[None, :])
                   + diag[:, None] + diag[None, :] + 2 * G)
            err[i1, :] = err[:, i1] = 1e18
            err[i2, :] = err[:, i2] = 1e18
            a, b = np.unravel_index(np.argmin(err), err.shape)
        # true residual ||obs - (sig_a+sig_b)||^2 (the dropped const restored)
        resid = float(obs @ obs - 2 * (oh[a] + oh[b])
                      + diag[a] + diag[b] + 2 * G[a, b])
        resids.append((s, max(resid, 0.0)))
        k_est = int(ocop[a] + ocop[b])
        by_k[k_true] += 1
        if k_est == k_true:
            correct += 1
            by_k_ok[k_true] += 1
        else:
            misses.append((s, k_true, k_est))

    n = len(samp)
    print(f"\noverall: {correct}/{n} = {100.0 * correct / n:.1f}%")
    print(f"{'k_true':>7} {'n':>5} {'correct':>8} {'acc':>7}")
    for k in sorted(by_k):
        print(f"{k:7d} {by_k[k]:5d} {by_k_ok[k]:8d} "
              f"{100.0 * by_k_ok[k] / by_k[k]:6.1f}%")
    if misses:
        print(f"\nmisses ({len(misses)}):")
        for s, kt, ke in misses[:30]:
            print(f"  {s:10s} true {kt}  est {ke}")

    rv = sorted(r for _, r in resids)
    nb = len(basis)
    exact0 = sum(1 for r in rv if r == 0.0)
    print(f"\nmatch residual ||obs-(sig_a+sig_b)||^2 over {nb} nodes:")
    print(f"  exact (residual 0): {exact0}/{len(rv)} samples")
    for q, lbl in [(0, "min"), (len(rv)//2, "median"), (len(rv)-1, "max")]:
        r = rv[q]
        print(f"  {lbl:6s} resid={r:10.1f}  per-node RMS={ (r/nb)**0.5:.4f}")
    worst = sorted(resids, key=lambda x: -x[1])[:5]
    print("  largest-residual samples: "
          + ", ".join(f"{s}({r:.0f})" for s, r in worst))


if __name__ == "__main__":
    main()
