#!/usr/bin/env python3
"""Split a CLR read FASTQ into one FASTQ per resolved orange paralog copy.

Each barcode in the cluster prediction is a read id. Reads assigned to a real
copy (pred_cluster `wh_<n>`) go to `<out-dir>/<cluster>.fq`; reads with no
usable markers (`wh_NA`) and reads absent from the prediction are dropped. The
number of output pools is data-dependent (at most the clustering ploidy).

Inputs:
  --pred   pred.tsv from cluster_whatshap.py (col 0 barcode, col 4 pred_cluster)
  --reads  the CLR read FASTQ that was fed to the pipeline

Output: one FASTQ per copy under --out-dir, named <pred_cluster>.fq."""

import argparse
import os
import sys


def load_pred(path):
    """read id -> pred_cluster, dropping the wh_NA non-cluster."""
    assign = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5 or c[4] == "wh_NA":
                continue
            assign[c[0]] = c[4]
    return assign


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True,
                    help="pred.tsv from cluster_whatshap.py")
    ap.add_argument("--reads", required=True, help="CLR read FASTQ")
    ap.add_argument("--out-dir", required=True,
                    help="directory for the per-copy FASTQ pools")
    args = ap.parse_args()

    assign = load_pred(args.pred)
    os.makedirs(args.out_dir, exist_ok=True)

    handles, n_written = {}, {}
    n_reads = n_kept = 0
    with open(args.reads) as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            rec = head + fh.readline() + fh.readline() + fh.readline()
            n_reads += 1
            rid = head[1:].split()[0]
            cl = assign.get(rid)
            if cl is None:
                continue
            if cl not in handles:
                handles[cl] = open(os.path.join(args.out_dir, f"{cl}.fq"), "w")
                n_written[cl] = 0
            handles[cl].write(rec)
            n_written[cl] += 1
            n_kept += 1
    for h in handles.values():
        h.close()

    for cl in sorted(n_written):
        sys.stderr.write(f"[INFO] split_pools: {cl} {n_written[cl]} reads\n")
    unassigned = n_reads - n_kept
    if unassigned:
        sys.stderr.write(f"[WARN] split_pools: {unassigned} reads not in any "
                         f"copy (wh_NA or absent from pred)\n")
    sys.stderr.write(f"[OK] split_pools: wrote {len(handles)} pools "
                     f"({n_kept}/{n_reads} reads) to {args.out_dir}\n")


if __name__ == "__main__":
    main()
