#!/usr/bin/env python3
"""Score resolved orange copies against the truth haplotypes.

Aligns every copy contig (minimap2, asm20) to both truth haplotype FASTAs and
keeps, per copy, the haplotype it matches best. Surfaces inflated copies
(length ratio >> 1, e.g. cyclic over-traversal) and partial / divergent ones.

Truth haplotypes are the region assemblies the reads were simulated from.

Output schema (--out): copy_id, truth_hap, copy_bp, truth_span_bp, len_ratio,
identity_pct, aln_pct, strand, note."""

import argparse
import os
import subprocess
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


def minimap(minimap2, target, query_path, preset="asm20"):
    """Run minimap2, return {query name -> best parsed PAF row}."""
    out = subprocess.run(
        [minimap2, "-c", "-x", preset, "--secondary=no", target, query_path],
        capture_output=True, text=True, check=True).stdout
    best = {}
    for ln in out.splitlines():
        f = ln.split("\t")
        if len(f) < 12:
            continue
        q = f[0]
        row = dict(qname=q, qlen=int(f[1]), qs=int(f[2]), qe=int(f[3]),
                   strand=f[4], ts=int(f[7]), te=int(f[8]),
                   nmatch=int(f[9]), alnlen=int(f[10]))
        row["div"] = next((float(t[5:]) for t in f[12:]
                           if t.startswith("de:f:")), None)
        if q not in best or row["nmatch"] > best[q]["nmatch"]:
            best[q] = row
    return best


def identity(row):
    if row["div"] is not None:
        return (1 - row["div"]) * 100
    return 100.0 * row["nmatch"] / max(row["alnlen"], 1)


def classify(ratio, alnf):
    if ratio > 1.3:
        return "INFLATED"
    if ratio < 0.7 or alnf < 70:
        return "PARTIAL"
    return "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--copies", required=True,
                    help="copies.fasta from collect_copies.py")
    ap.add_argument("--truth-h1", required=True,
                    help="haplotype-1 truth region assembly FASTA")
    ap.add_argument("--truth-h2", required=True,
                    help="haplotype-2 truth region assembly FASTA")
    ap.add_argument("--minimap2", required=True, help="minimap2 binary path")
    ap.add_argument("--out", required=True, help="per-copy accuracy TSV")
    args = ap.parse_args()

    copies = read_fasta(args.copies)
    aln = {"h1": minimap(args.minimap2, args.truth_h1, args.copies),
           "h2": minimap(args.minimap2, args.truth_h2, args.copies)}

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n_ok = n_bad = n_unaligned = 0
    with open(args.out, "w") as out:
        out.write("copy_id\ttruth_hap\tcopy_bp\ttruth_span_bp\tlen_ratio\t"
                  "identity_pct\taln_pct\tstrand\tnote\n")
        for cid in sorted(copies):
            seq = copies[cid]
            r1, r2 = aln["h1"].get(cid), aln["h2"].get(cid)
            cand = [(h, r) for h, r in (("h1", r1), ("h2", r2)) if r]
            if not cand:
                out.write(f"{cid}\t.\t{len(seq)}\t.\t.\t.\t.\t.\tUNALIGNED\n")
                n_unaligned += 1
                continue
            hap, r = max(cand, key=lambda hr: hr[1]["nmatch"])
            tspan = r["te"] - r["ts"]
            ratio = len(seq) / max(tspan, 1)
            alnf = 100.0 * (r["qe"] - r["qs"]) / max(len(seq), 1)
            note = classify(ratio, alnf)
            n_ok += note == "ok"
            n_bad += note != "ok"
            out.write(f"{cid}\t{hap}\t{len(seq)}\t{tspan}\t{ratio:.3f}\t"
                      f"{identity(r):.2f}\t{alnf:.1f}\t{r['strand']}\t{note}\n")

    sys.stderr.write(f"[OK] compare_to_truth.py: {len(copies)} copies "
                     f"({n_ok} ok, {n_bad} inflated/partial, "
                     f"{n_unaligned} unaligned) -> {args.out}\n")


if __name__ == "__main__":
    main()
