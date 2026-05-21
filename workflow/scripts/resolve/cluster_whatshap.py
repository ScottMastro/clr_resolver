#!/usr/bin/env python3
"""Deconvolve orange barcodes with WhatsHap polyphase's actual engine.

Adapter: our per-barcode marker matrix maps 1:1 onto polyphase's read model
(a Read = name + [(position, allele, quality)]). Each barcode becomes one
Read, each informative bubble a variant position. We then call polyphase's
own Phase I (cluster editing) and, optionally, Phase II (haplotype threading).

Run with the whatshap-env interpreter:
  /home/scott/anaconda3/envs/whatshap-env/bin/python cluster_whatshap.py ...
"""

import argparse
import sys
from collections import defaultdict

from whatshap.core import Read, ReadSet
from whatshap.polyphase.solver import (
    scoreReadset, ClusterEditingSolver, AlleleMatrix,
)
from whatshap.polyphase.threading import run_threading


def load_matrix(path):
    mat, src = {}, {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            bc = c[0]
            src[bc] = c[1]
            row = {i: int(v) for i, v in enumerate(c[4:]) if v in ("0", "1")}
            mat[bc] = row
    return mat, src


def build_readset(mat):
    """One Read per barcode; marker index = variant position, allele 0/1."""
    rs = ReadSet()
    order = []
    for bc in sorted(mat):
        row = mat[bc]
        if not row:
            continue
        positions = sorted(row)
        r = Read(bc, 60, 0, 0, positions[0], bc)   # name, mapq, src, sample, refstart, BX
        for p in positions:
            r.add_variant(p, row[p], 30)
        r.sort()
        rs.add(r)
        order.append(bc)
    rs.sort()
    return rs, order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ploidy", type=int, default=6)
    ap.add_argument("--min-overlap", type=int, default=2)
    ap.add_argument("--err", type=float, default=0.07)
    ap.add_argument("--bundle-edges", action="store_true")
    ap.add_argument("--thread", action="store_true",
                    help="run Phase II threading to collapse clusters to "
                         "--ploidy haplotypes")
    ap.add_argument("--collapse", action="store_true",
                    help="collapse Phase I clusters to --ploidy copies via "
                         "centroid agreement instead of threading")
    ap.add_argument("--refine", type=int, default=30,
                    help="EM iterations for --collapse")
    args = ap.parse_args()

    mat, src = load_matrix(args.matrix)
    rs, order = build_readset(mat)
    # ReadSet.sort() may reorder; recover barcode order from read names.
    order = [rs[i].name for i in range(len(rs))]
    print(f"[readset] {len(rs)} barcodes as reads", file=sys.stderr)

    am = AlleleMatrix(rs)
    n_pos = len(am.getPositions())
    print(f"[matrix] {n_pos} variant positions", file=sys.stderr)

    # Phase I -- cluster editing
    sim = scoreReadset(am, args.min_overlap, args.ploidy, args.err)
    solver = ClusterEditingSolver(sim, args.bundle_edges)
    clustering = solver.run()
    for i in range(sum(len(c) for c in clustering), len(rs)):
        clustering.append([i])
    sizes = sorted((len(c) for c in clustering), reverse=True)
    print(f"[cluster-editing] {len(clustering)} clusters; "
          f"sizes={sizes[:12]}", file=sys.stderr)

    if args.collapse:
        # polyphase Phase I clusters are excellent but over-segmented. Collapse
        # like k-modes: the `ploidy` largest clusters seed full centroids,
        # then every barcode is assigned to the centroid it best agrees with;
        # a few EM iterations stitch the fragments back together.
        clustering.sort(key=len, reverse=True)
        bc_of = {ri: order[ri] for ri in range(len(rs))}

        def full_consensus(members):
            tally = defaultdict(lambda: [0, 0])
            for ri in members:
                for m, v in mat[bc_of[ri]].items():
                    tally[m][v] += 1
            return {m: (1 if t[1] > t[0] else 0)
                    for m, t in tally.items() if t[0] + t[1] > 0}

        def agree(row, cent):
            n = ok = 0
            for m, v in row.items():
                cv = cent.get(m)
                if cv is None:
                    continue
                n += 1
                ok += (v == cv)
            return (ok / n, n) if n else (0.0, 0)

        cents = [full_consensus(c) for c in clustering[:args.ploidy]]
        cur = {}
        for _ in range(args.refine + 1):
            new = {}
            for bc in sorted(mat):
                row = mat[bc]
                bh, ba, bn = 0, -1.0, -1
                for ci, ce in enumerate(cents):
                    a, n = agree(row, ce)
                    if a > ba or (a == ba and n > bn):
                        bh, ba, bn = ci, a, n
                new[bc] = bh
            if new == cur:
                break
            cur = new
            grp = defaultdict(list)
            bc2ri = {order[ri]: ri for ri in range(len(rs))}
            for bc, ci in cur.items():
                if bc in bc2ri:
                    grp[ci].append(bc2ri[bc])
            for ci in range(len(cents)):
                if grp[ci]:
                    cents[ci] = full_consensus(grp[ci])
        assign = {bc: f"wh_{ci}" for bc, ci in cur.items()}
    elif args.thread:
        # Phase II -- thread `ploidy` haplotypes through the clusters.
        # Genotype per position: allele-1 dosage estimated from the pileup,
        # rounded to the ploidy. distrust_genotypes lets threading override.
        gt = []
        for p in range(n_pos):
            c0 = c1 = 0
            for ri in range(len(rs)):
                a = am.getAllele(ri, p)
                if a == 0:
                    c0 += 1
                elif a == 1:
                    c1 += 1
            tot = c0 + c1
            d = round(args.ploidy * c1 / tot) if tot else 0
            gt.append({0: args.ploidy - d, 1: d})
        threads, haps = run_threading(am, clustering, args.ploidy, gt,
                                      distrust_genotypes=True)
        # haps[h] = the threaded consensus allele sequence of haplotype h,
        # indexed by position. Assign every barcode to the haplotype whose
        # sequence it best agrees with over the positions it covers.
        positions = list(am.getPositions())
        pidx = {p: i for i, p in enumerate(positions)}
        assign = {}
        for bc in sorted(mat):
            row = mat[bc]
            best_h, best_a, best_n = -1, -1.0, -1
            for h in range(args.ploidy):
                seq = haps[h]
                ok = n = 0
                for mi, al in row.items():
                    hv = seq[pidx[mi]]
                    if hv < 0:
                        continue
                    n += 1
                    ok += (al == hv)
                a = ok / n if n else 0.0
                if a > best_a or (a == best_a and n > best_n):
                    best_h, best_a, best_n = h, a, n
            assign[bc] = f"wh_{best_h}" if best_h >= 0 else "wh_NA"
    else:
        assign = {}
        for ci, cl in enumerate(clustering):
            for ri in cl:
                assign[order[ri]] = f"wh_{ci}"

    with open(args.out, "w") as out:
        out.write("barcode\tsource_hap\tstruct_class\tn_markers\tpred_cluster\n")
        for bc in sorted(mat):
            out.write(f"{bc}\t{src[bc]}\tO\t{len(mat[bc])}\t"
                      f"{assign.get(bc, 'wh_NA')}\n")
    print(f"[OK] {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
