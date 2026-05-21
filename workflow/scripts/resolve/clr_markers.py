#!/usr/bin/env python3
"""CLR-aware bubble allele calling — walk the read, don't trust the route.

discover_markers.py calls a bubble allele from *which interior node the path
runs through*. For a CLR read that is the aligner's forced choice: when the
informative base is deleted (the common CLR error) or substituted to a base
that is neither allele, the read carries no real evidence and the aligner
routes it onto a branch arbitrarily — a ~9% per-call miscall.

Here we instead walk the read's CIGAR (GAF `cg:Z:` tag) against its path. A
bubble call is emitted only when the interior allele node aligns as an exact
match (`=`). If the node's base is a deletion (`D`) or mismatch (`X`) in the
read, the read had no evidence — we emit nothing (a missing call) rather than
a forced 0/1. Forced-random noise becomes honest missingness, which the
downstream phi / EM steps handle natively.

Outputs match discover_markers.py exactly, so build_matrix.py is unchanged:
  --out-markers  bubble_id, src, snk, allele0, allele1, len0, len1,
                 n0, n1, kind, informative
  --out-calls    read_id, mate, barcode, source_hap, bubble_id, allele
"""

import argparse
import json
import re
import sys
from collections import defaultdict

_STEP = re.compile(r"([<>])(\d+)")
_CIG = re.compile(r"(\d+)([=XIDM])")


def load_orange(path, color):
    """Nodes of the target colour. color='ALL' = every coloured node, so
    bubbles in any colour become markers (used to add cross-SD-breakpoint
    flanking context when untangling one colour's copies)."""
    s = set()
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if c[1] == color or (color == "ALL" and c[1] != "."):
                s.add(c[0])
    return s


def load_node_len(gfa):
    ln = {}
    with open(gfa) as fh:
        for line in fh:
            if not line.startswith("S\t"):
                continue
            p = line.rstrip("\n").split("\t")
            seq = p[2] if len(p) > 2 else "*"
            if seq != "*":
                ln[p[1]] = len(seq)
            else:
                m = re.search(r"\bLN:i:(\d+)\b", line)
                ln[p[1]] = int(m.group(1)) if m else 0
    return ln


def path_ops(cigar, pstart, path_len):
    """Op per path base over [pstart, pstart+consumed): '=', 'X' or 'D'.
    Returns a list indexed by absolute path coordinate (None outside)."""
    ops = [None] * path_len
    pos = pstart
    for n, op in _CIG.findall(cigar):
        n = int(n)
        if op == "I":                       # consumes read only
            continue
        sym = "=" if op == "M" else op       # =, X, D consume path
        for _ in range(n):
            if 0 <= pos < path_len:
                ops[pos] = sym
            pos += 1
    return ops


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gaf", required=True)
    ap.add_argument("--bubbles", required=True)
    ap.add_argument("--node-color", required=True)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--color", default="O")
    ap.add_argument("--out-markers", required=True)
    ap.add_argument("--out-calls", required=True)
    ap.add_argument("--min-allele-obs", type=int, default=3)
    args = ap.parse_args()

    orange = load_orange(args.node_color, args.color)
    node_len = load_node_len(args.gfa)
    chains = json.load(open(args.bubbles))

    bubbles, inside_idx = {}, {}
    for chain in chains.values():
        for b in chain.get("bubbles", []):
            if b.get("type") != "simple":
                continue
            ins = b.get("inside", [])
            if len(ins) != 2 or not all(n in orange for n in ins):
                continue
            bid = b["id"]
            ends = b.get("ends", ["", ""])
            l0, l1 = node_len.get(ins[0], 0), node_len.get(ins[1], 0)
            kind = "snp" if l0 == 1 and l1 == 1 else (
                "indel" if l0 != l1 else "mnp")
            bubbles[bid] = dict(src=ends[0], snk=ends[1], a0=ins[0],
                                a1=ins[1], l0=l0, l1=l1, kind=kind, n0=0, n1=0)
            inside_idx[ins[0]] = (bid, 0)
            inside_idx[ins[1]] = (bid, 1)
    print(f"[info] {len(bubbles)} simple biallelic bubbles inside color "
          f"{args.color}", file=sys.stderr)

    # best GAF alignment per read (longest block, col 11 / idx 10)
    best = {}
    with open(args.gaf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            rid, blocklen = f[0], int(f[10])
            if rid not in best or blocklen > best[rid][0]:
                cg = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), "")
                best[rid] = (blocklen, f[5], int(f[7]), int(f[8]), cg)

    calls = []
    n_conf = n_forced = n_oob = 0
    for rid, (_, pathstr, pstart, pend, cigar) in best.items():
        if not cigar:
            continue
        steps = [(o, n) for o, n in _STEP.findall(pathstr)]
        # cumulative path offset of each step
        offs, cum = [], 0
        for _, n in steps:
            offs.append(cum)
            cum += node_len.get(n, 0)
        ops = path_ops(cigar, pstart, cum)

        hit = {}        # bubble_id -> set of (allele, confident)
        for (o, n), off in zip(steps, offs):
            if n not in inside_idx:
                continue
            bid, al = inside_idx[n]
            nlen = node_len.get(n, 0)
            seg = ops[off:off + nlen]
            if nlen == 0 or any(s is None for s in seg):
                conf = None                   # outside aligned span
            else:
                conf = all(s == "=" for s in seg)
            hit.setdefault(bid, set()).add((al, conf))

        src = rid.split("_", 1)[0]
        for bid, obs in hit.items():
            alleles = {al for al, _ in obs}
            if len(alleles) != 1:             # walked both branches
                continue
            al = next(iter(alleles))
            conf = any(c for _, c in obs)
            oob = all(c is None for _, c in obs)
            if conf:
                n_conf += 1
                calls.append((rid, "1", rid, src, bid, al))
                bubbles[bid]["n0" if al == 0 else "n1"] += 1
            elif oob:
                n_oob += 1                    # node outside aligned span
            else:
                n_forced += 1                 # D/X at the informative base

    with open(args.out_markers, "w") as out:
        out.write("bubble_id\tsrc\tsnk\tallele0\tallele1\tlen0\tlen1\t"
                  "n0\tn1\tkind\tinformative\n")
        n_inf = 0
        for bid, b in sorted(bubbles.items()):
            inf = (b["n0"] >= args.min_allele_obs and
                   b["n1"] >= args.min_allele_obs)
            n_inf += inf
            out.write(f"{bid}\t{b['src']}\t{b['snk']}\t{b['a0']}\t{b['a1']}\t"
                      f"{b['l0']}\t{b['l1']}\t{b['n0']}\t{b['n1']}\t"
                      f"{b['kind']}\t{int(inf)}\n")

    with open(args.out_calls, "w") as out:
        out.write("read_id\tmate\tbarcode\tsource_hap\tbubble_id\tallele\n")
        for r in calls:
            out.write("\t".join(map(str, r)) + "\n")

    total = n_conf + n_forced
    print(f"[OK] {args.out_markers}: {n_inf} informative", file=sys.stderr)
    print(f"[OK] {args.out_calls}: {n_conf} confident calls kept; "
          f"{n_forced} forced (D/X at the base) dropped to missing "
          f"({n_forced/max(total,1):.1%}); {n_oob} out-of-span",
          file=sys.stderr)


if __name__ == "__main__":
    main()
