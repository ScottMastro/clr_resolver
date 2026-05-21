#!/usr/bin/env python3
"""Combine one benchmark sample's metrics into a summary row and a report.

Reports the three benchmark questions:
  (a) read phasing accuracy  -- from score.py
  (b) copy count correct     -- resolved copies vs distinct truth copies
  (c) bp accuracy            -- mean per-copy identity and length ratio

Inputs:
  --phasing          phasing.tsv from score.py (metric/value)
  --bp               bp_accuracy.tsv from compare_to_truth.py
  --copies-manifest  copies.tsv from collect_copies.py
  --truth            truth.tsv from clr_truth.py

Output schema (--out-tsv): sample, phasing_accuracy, reads_scored,
copies_found, copies_expected, copy_count_correct, bp_identity_mean,
bp_len_ratio_mean, bp_copies_ok."""

import argparse
import os
import sys


def load_kv(path):
    """metric/value TSV -> dict."""
    d = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2:
                d[c[0]] = c[1]
    return d


def load_rows(path):
    """TSV with header -> list of dicts."""
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]


def mean(vals):
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def as_float(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def fmt(x, nd=3):
    return f"{x:.{nd}f}" if x is not None else "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phasing", required=True, help="score.py metrics TSV")
    ap.add_argument("--bp", required=True,
                    help="compare_to_truth.py per-copy TSV")
    ap.add_argument("--copies-manifest", required=True,
                    help="collect_copies.py copies.tsv")
    ap.add_argument("--truth", required=True, help="clr_truth.py truth.tsv")
    ap.add_argument("--sample", required=True, help="benchmark sample id")
    ap.add_argument("--out-tsv", required=True, help="one-row summary TSV")
    ap.add_argument("--out-txt", required=True, help="human-readable summary")
    args = ap.parse_args()

    phasing = load_kv(args.phasing)
    acc = as_float(phasing.get("accuracy"))
    reads_scored = phasing.get("reads_scored", "NA")

    manifest = load_rows(args.copies_manifest)
    copies_found = sum(1 for r in manifest
                       if as_float(r.get("contig_len")) not in (None, 0.0))

    truth_copies = set()
    with open(args.truth) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 5 and c[4] != ".":
                truth_copies.add(c[4])
    copies_expected = len(truth_copies)
    copy_count_correct = copies_found == copies_expected

    bp = load_rows(args.bp)
    ident = mean([as_float(r.get("identity_pct")) for r in bp])
    lenr = mean([as_float(r.get("len_ratio")) for r in bp])
    bp_ok = sum(1 for r in bp if r.get("note") == "ok")

    os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)
    with open(args.out_tsv, "w") as out:
        out.write("sample\tphasing_accuracy\treads_scored\tcopies_found\t"
                  "copies_expected\tcopy_count_correct\tbp_identity_mean\t"
                  "bp_len_ratio_mean\tbp_copies_ok\n")
        out.write(f"{args.sample}\t{fmt(acc, 4)}\t{reads_scored}\t"
                  f"{copies_found}\t{copies_expected}\t"
                  f"{int(copy_count_correct)}\t{fmt(ident, 2)}\t"
                  f"{fmt(lenr, 3)}\t{bp_ok}\n")

    os.makedirs(os.path.dirname(args.out_txt) or ".", exist_ok=True)
    with open(args.out_txt, "w") as out:
        out.write(f"clr_resolve benchmark -- {args.sample}\n")
        out.write("=" * 48 + "\n\n")
        out.write(f"(a) read phasing accuracy : {fmt(acc, 3)}  "
                  f"({reads_scored} reads scored)\n")
        verdict = "CORRECT" if copy_count_correct else "WRONG"
        out.write(f"(b) copy count            : {copies_found} resolved / "
                  f"{copies_expected} expected  [{verdict}]\n")
        out.write(f"(c) bp accuracy           : {fmt(ident, 2)}% mean identity"
                  f", {fmt(lenr, 3)} mean length ratio\n")
        out.write(f"                            {bp_ok}/{len(bp)} copies "
                  f"within expected length bounds\n")

    sys.stderr.write(f"[OK] summarize_benchmark: {args.sample} "
                     f"acc={fmt(acc, 3)} copies={copies_found}/"
                     f"{copies_expected} ident={fmt(ident, 2)}\n")


if __name__ == "__main__":
    main()
