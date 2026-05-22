#!/usr/bin/env python3
"""Evaluate orange copy-number estimation on real read GAFs (predictor A).

For each sim_<sample> GAF, build the observed read-depth vector over the
anchor-node basis, drop the sample's own haplotypes from the panel (LOO), and
match to the best remaining haplotype pair. k = matched pair's painted copy
count. Reports cosine and least-squares matches, all- and orange-basis,
against the painting-truth k.

Output: per-sample table + accuracy to stdout."""

import argparse
import glob
import os
import re
import sys
from collections import Counter

import numpy as np

_SEG = re.compile(r"(\d+)[+-]")
_STEP = re.compile(r"[<>](\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--node-color", required=True)
    ap.add_argument("--gaf-dir", required=True, help="dir of sim_<sample>.gaf")
    ap.add_argument("--max-mult", type=int, default=8)
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
    in_all = [n for n in nlen if sum(n in cnt[nm] for nm in names) == npath]
    maxmult = {n: max(cnt[nm].get(n, 0) for nm in names) for n in in_all}
    anchors = [n for n in in_all if maxmult[n] <= args.max_mult]
    sc = [n for n in anchors if maxmult[n] == 1]          # depth-scale nodes
    obasis = [n for n in anchors if col.get(n, ".") == "O"]
    # full orange set: every orange node (incl. copy-variable, non-anchor),
    # VNTR-capped -- carries the copy-distinguishing signal anchors discard
    omax = {}
    for nm in names:
        for n, m in cnt[nm].items():
            if col.get(n, ".") == "O" and m > omax.get(n, 0):
                omax[n] = m
    oall = [n for n in omax if omax[n] <= args.max_mult]
    sys.stderr.write(f"[INFO] {npath} haplotypes; anchors {len(anchors)} "
                     f"(orange-anchor {len(obasis)}, orange-all {len(oall)}); "
                     f"single-copy {len(sc)}\n")

    def sigmat(basis):
        idx = {n: i for i, n in enumerate(basis)}
        M = np.zeros((npath, len(basis)))
        for i, nm in enumerate(names):
            for n, m in cnt[nm].items():
                if n in idx:
                    M[i, idx[n]] = m
        return M, idx

    SIG_all, idx_all = sigmat(anchors)
    SIG_org, idx_org = sigmat(obasis)
    SIG_oa, idx_oa = sigmat(oall)
    ocop = np.array([ocopies(paths[nm]) for nm in names])
    rowof = {nm.split("#")[0]: [] for nm in names}
    for i, nm in enumerate(names):
        rowof[nm.split("#")[0]].append(i)

    def gaf_depth(path):
        d = {}
        for ln in open(path):
            f = ln.split("\t")
            if len(f) < 6:
                continue
            for n in _STEP.findall(f[5]):
                d[n] = d.get(n, 0) + 1
        return d

    def match(SIG, idx, depth, samp, metric, d):
        obs = np.zeros(SIG.shape[1])
        for n, i in idx.items():
            obs[i] = depth.get(n, 0)
        G = SIG @ SIG.T
        diag = np.diag(G)
        rows = rowof.get(samp, [])
        if metric == "cosine":
            oh = SIG @ obs
            den = np.sqrt(np.maximum(diag[:, None] + diag[None, :] + 2 * G,
                                     1e-9)) * np.sqrt(max(obs @ obs, 1e-9))
            sc_ = (oh[:, None] + oh[None, :]) / den
            for r in rows:
                sc_[r, :] = sc_[:, r] = -1
            a, b = np.unravel_index(np.argmax(sc_), sc_.shape)
        else:
            t = obs / d
            oh = SIG @ t
            err = (-2 * (oh[:, None] + oh[None, :])
                   + diag[:, None] + diag[None, :] + 2 * G)
            for r in rows:
                err[r, :] = err[:, r] = 1e18
            a, b = np.unravel_index(np.argmin(err), err.shape)
        return int(ocop[a] + ocop[b])

    gafs = sorted(glob.glob(os.path.join(args.gaf_dir, "sim_*.gaf")))
    print(f"{'sample':12s} {'truek':>5s} {'oCOS':>5s} {'oLSQ':>5s} "
          f"{'faCOS':>6s} {'faLSQ':>6s}")
    tally = Counter()
    for g in gafs:
        samp = os.path.basename(g)[4:-4]            # sim_X.gaf -> X
        if samp not in rowof:
            continue
        rows = rowof[samp]
        truek = int(sum(ocop[r] for r in rows))
        dep = gaf_depth(g)
        d = np.median([dep.get(n, 0) for n in sc]) / 2 or 1.0
        oc_ = match(SIG_org, idx_org, dep, samp, "cosine", d)
        ol = match(SIG_org, idx_org, dep, samp, "lsq", d)
        fac = match(SIG_oa, idx_oa, dep, samp, "cosine", d)
        fal = match(SIG_oa, idx_oa, dep, samp, "lsq", d)
        for tag, est in (("oCOS", oc_), ("oLSQ", ol),
                         ("faCOS", fac), ("faLSQ", fal)):
            tally[tag] += (est == truek)
        print(f"{samp:12s} {truek:5d} {oc_:5d} {ol:5d} {fac:6d} {fal:6d}")
    n = sum(1 for g in gafs if os.path.basename(g)[4:-4] in rowof)
    print(f"\naccuracy over {n} read sets:")
    for tag in ("oCOS", "oLSQ", "faCOS", "faLSQ"):
        print(f"  {tag}: {tally[tag]}/{n} = {100.0*tally[tag]/n:.1f}%")


if __name__ == "__main__":
    main()
