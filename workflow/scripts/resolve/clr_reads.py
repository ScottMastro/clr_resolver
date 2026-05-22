#!/usr/bin/env python3
"""Adapt a CLR graph-GAF into the deconv pipeline's per-read reads table.

The deconv copy-axis pipeline (extract_reads -> discover_markers ->
build_matrix -> cluster_sda) treats one *molecule* per row, keyed by a 10x
barcode. CLR has no barcode: each long read IS the molecule — it spans many
PSVs directly, which is the setting SDA was originally built for. So we emit
one row per CLR read with barcode == read_id.

Input GAF read names are `h1_S1_<n>` / `h2_S1_<n>`; the `h1`/`h2` prefix is
the simulation truth haplotype. A read may have several GAF alignments — we
keep the one with the longest alignment block (col 11).

Output (TSV, extract_reads.py-compatible 12 columns):
  read_id, mate, source_hap, hap_idx, mol_start, mol_end,
  read_start, read_end, barcode, n_nodes, n_C_nodes, path
mol_*/read_* are unused downstream for CLR and emitted as 0; barcode == the
read id (one molecule per read); source_hap is the h1/h2 truth tag.
"""

import argparse
import gzip
import re
import sys

_STEP = re.compile(r"\d+")


def open_text(path):
    """Open a text file, transparently decompressing a .gz by extension."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_color(path, target):
    in_c = set()
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[1] == target:
                in_c.add(c[0])
    return in_c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaf", required=True)
    ap.add_argument("--node-color", required=True)
    ap.add_argument("--color", default="O")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-c-nodes", type=int, default=2,
                    help="keep reads touching >= this many color-C nodes")
    args = ap.parse_args()

    in_c = load_color(args.node_color, args.color)
    print(f"[info] color {args.color}: {len(in_c)} nodes", file=sys.stderr)

    # best alignment per read = longest alignment block (GAF col 11, idx 10)
    best = {}
    n_aln = 0
    with open_text(args.gaf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            n_aln += 1
            rid, path, blocklen = f[0], f[5], int(f[10])
            if rid not in best or blocklen > best[rid][0]:
                best[rid] = (blocklen, path)

    n_kept = 0
    with open(args.out, "w") as out:
        out.write("read_id\tmate\tsource_hap\thap_idx\tmol_start\tmol_end\t"
                  "read_start\tread_end\tbarcode\tn_nodes\tn_C_nodes\tpath\n")
        for rid, (_, path) in sorted(best.items()):
            nodes = _STEP.findall(path)
            n_in_c = sum(1 for n in nodes if n in in_c)
            if n_in_c < args.min_c_nodes:
                continue
            n_kept += 1
            src = rid.split("_", 1)[0]            # h1 / h2 truth tag
            out.write(f"{rid}\t1\t{src}\t0\t0\t0\t0\t0\t{rid}\t"
                      f"{len(nodes)}\t{n_in_c}\t{','.join(nodes)}\n")

    print(f"[OK] {args.out}: {n_aln} alignments -> {len(best)} reads, "
          f"{n_kept} in color {args.color}", file=sys.stderr)


if __name__ == "__main__":
    main()
