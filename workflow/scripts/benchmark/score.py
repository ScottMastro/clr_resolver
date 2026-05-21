#!/usr/bin/env python3
"""Score a deconvolution against the truth set.

Predicted clusters are unlabeled, so we find the best one-to-one matching
between predicted clusters and truth copies (greedy on the contingency
table), then report:

  - accuracy   = fraction of barcodes whose predicted cluster maps to their
                 true copy
  - purity     = per predicted cluster, the fraction from its dominant copy
  - recovery   = per truth copy, fraction of its barcodes captured by the
                 matched cluster

The per-cluster table is printed to stdout; --out-tsv writes the headline
metrics as a key/value table.

Output schema (--out-tsv): metric, value -- rows reads_scored, pred_clusters,
truth_copies, accuracy."""

import argparse
import os
import sys
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True,
                    help="pred.tsv (barcode -> pred_cluster)")
    ap.add_argument("--truth", required=True,
                    help="truth.tsv (barcode -> truth_copy)")
    ap.add_argument("--out-tsv", required=True,
                    help="headline metrics as a key/value TSV")
    args = ap.parse_args()

    truth = {}
    with open(args.truth) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            truth[c[0]] = c[4]                # barcode -> truth_copy

    pred = {}
    with open(args.pred) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            pred[c[0]] = c[4]                 # barcode -> pred_cluster

    bcs = [b for b in pred if b in truth]
    # contingency: pred_cluster x truth_copy
    cont = defaultdict(lambda: defaultdict(int))
    for b in bcs:
        cont[pred[b]][truth[b]] += 1

    # greedy best matching pred->truth
    pairs = []
    for pc, row in cont.items():
        for tc, n in row.items():
            pairs.append((n, pc, tc))
    pairs.sort(reverse=True)
    pc2tc, used_tc = {}, set()
    for n, pc, tc in pairs:
        if pc in pc2tc or tc in used_tc:
            continue
        pc2tc[pc] = tc
        used_tc.add(tc)

    correct = sum(1 for b in bcs if pc2tc.get(pred[b]) == truth[b])
    acc = correct / len(bcs) if bcs else 0.0
    n_truth_copies = len(set(truth[b] for b in bcs))

    print(f"barcodes scored: {len(bcs)}")
    print(f"predicted clusters: {len(cont)}   truth copies: {n_truth_copies}")
    print(f"ACCURACY (correct copy): {correct}/{len(bcs)} = {acc:.3f}")
    print()
    print(f"{'pred_cluster':16} {'-> truth':14} {'size':>5} {'purity':>7} "
          f"{'recovery':>9}")
    tc_total = defaultdict(int)
    for b in bcs:
        tc_total[truth[b]] += 1
    for pc in sorted(cont, key=lambda p: -sum(cont[p].values())):
        row = cont[pc]
        size = sum(row.values())
        tc = pc2tc.get(pc, "-")
        dom = row.get(tc, 0) if tc != "-" else max(row.values())
        purity = dom / size if size else 0.0
        recov = dom / tc_total[tc] if tc != "-" and tc_total[tc] else 0.0
        comp = ",".join(f"{t}:{n}" for t, n in
                        sorted(row.items(), key=lambda x: -x[1]))
        print(f"{pc:16} {tc:14} {size:5} {purity:7.3f} {recov:9.3f}  {comp}")

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    with open(args.out_tsv, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"reads_scored\t{len(bcs)}\n")
        out.write(f"pred_clusters\t{len(cont)}\n")
        out.write(f"truth_copies\t{n_truth_copies}\n")
        out.write(f"accuracy\t{acc:.4f}\n")
    sys.stderr.write(f"[OK] score.py: wrote {args.out_tsv} (accuracy {acc:.3f})\n")


if __name__ == "__main__":
    main()
