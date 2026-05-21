#!/usr/bin/env python3
"""Gather per-pool Flye assemblies into one multi-FASTA of orange copies.

A Flye pool normally yields one contig; it may occasionally yield several
(mild fragmentation) or none (assembly failure). Every contig is emitted and
recorded in the manifest, so a single bad pool does not sink the dataset.

Inputs: one or more Flye assembly.fasta paths, each laid out as
  .../flye/<dataset>/<pool>/assembly.fasta
so the pool id is the parent directory name.

Output schema:
  --out-fasta     copy contigs, headers `<pool>.<n>`
  --out-manifest  TSV: copy_id, source_pool, contig_len, n_contigs_in_pool"""

import argparse
import os
import sys


def read_fasta(path):
    """header -> sequence."""
    seqs, name = {}, None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                name = ln[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assemblies", required=True, nargs="+",
                    help="Flye assembly.fasta paths, one per pool")
    ap.add_argument("--out-fasta", required=True,
                    help="merged copy FASTA")
    ap.add_argument("--out-manifest", required=True,
                    help="copy manifest TSV")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_fasta) or ".", exist_ok=True)
    rows, n_copies = [], 0
    with open(args.out_fasta, "w") as out:
        for path in sorted(args.assemblies):
            pool = os.path.basename(os.path.dirname(path))
            seqs = read_fasta(path) if os.path.exists(path) else {}
            if not seqs:
                sys.stderr.write(f"[WARN] collect_copies: pool {pool} produced "
                                 f"no contig\n")
                rows.append((f"{pool}.0", pool, 0, 0))
                continue
            for i, (_, seq) in enumerate(sorted(seqs.items())):
                cid = f"{pool}.{i}"
                out.write(f">{cid}\n{seq}\n")
                rows.append((cid, pool, len(seq), len(seqs)))
                n_copies += 1

    os.makedirs(os.path.dirname(args.out_manifest) or ".", exist_ok=True)
    with open(args.out_manifest, "w") as mf:
        mf.write("copy_id\tsource_pool\tcontig_len\tn_contigs_in_pool\n")
        for cid, pool, ln, n in rows:
            mf.write(f"{cid}\t{pool}\t{ln}\t{n}\n")

    sys.stderr.write(f"[OK] collect_copies: wrote {n_copies} copies from "
                     f"{len(args.assemblies)} pools to {args.out_fasta}\n")


if __name__ == "__main__":
    main()
