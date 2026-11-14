#!/usr/bin/env python3
"""Periodicity gate: a PASS / WARN / FAIL verdict from the output of `ribotish quality`.

Inputs
  qual.txt  line 1: {read_length: count}; line 4: {read_length: [frame0, frame1, frame2]} over annotated CDSs
  para.py   `offdict = {length: P-site offset, ..., 'm0': {...}}`, the read lengths Ribo-TISH accepted as periodic

Rules
  FAIL  no accepted read length (no 3-nt periodicity)
  FAIL  accepted lengths hold less than --min-pct of the reads (default 50; calibrated on RNase I footprints, so
        lower it per dataset, with a reason, for nucleases that give broader length distributions)
  WARN  fewer than --min-reads reads in accepted lengths (default 5,000,000)
  WARN  an offset outside the plausible range for its length (10 to 14 nt for 25 to 33 nt reads, 13 to 18 nt for
        34 to 42 nt reads, 9 to 20 nt otherwise)

Exit code 1 on FAIL, 0 otherwise. Both input files are parsed as literals; nothing in them is executed.
"""
import argparse
import ast
import sys
from collections import Counter


def read_qual(path):
    with open(path) as fh:
        lines = [ln.strip() for ln in fh]
    if len(lines) < 4:
        sys.exit(f"{path}: expected at least 4 lines of ribotish quality output, found {len(lines)}")
    return ast.literal_eval(lines[0]), ast.literal_eval(lines[3])


def read_offsets(path):
    with open(path) as fh:
        tree = ast.parse(fh.read())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "offdict" for t in node.targets):
            return ast.literal_eval(node.value)
    return {}


def plausible(length, offset):
    if 25 <= length <= 33:
        return 10 <= offset <= 14
    if 34 <= length <= 42:
        return 13 <= offset <= 18
    return 9 <= offset <= 20


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("qual")
    ap.add_argument("para")
    ap.add_argument("--label", required=True)
    ap.add_argument("--min-reads", type=int, default=5_000_000)
    ap.add_argument("--min-pct", type=float, default=50.0)
    ap.add_argument("--tsv", help="write the verdict as a one-row table")
    args = ap.parse_args()

    counts, frames = read_qual(args.qual)
    accepted = {k: v for k, v in read_offsets(args.para).items() if k != "m0"}
    total = sum(counts.values())
    acc_reads = sum(counts.get(length, 0) for length in accepted)
    pct = 100.0 * acc_reads / total if total else 0.0

    print(f"periodicity gate: {args.label}")
    print(f"  reads in the tested lengths: {total:,}")
    print(f"  {'len':>4} {'reads':>12} {'%':>6} {'frame0':>10} {'frame1':>10} {'frame2':>10} {'best':>5} offset")
    for length in sorted(counts, key=lambda x: -counts[x]):
        f = frames.get(length, [0, 0, 0])
        best = max(f) / (sum(f) or 1)
        off = accepted.get(length)
        flag = "" if off is None or plausible(length, off) else "  implausible offset"
        print(f"  {length:>4} {counts[length]:>12,} {100 * counts[length] / (total or 1):>5.1f}% {f[0]:>10,} "
              f"{f[1]:>10,} {f[2]:>10,} {best:>5.2f} {'-' if off is None else off}{flag}")
    print(f"  accepted lengths and offsets: { {k: accepted[k] for k in sorted(accepted)} if accepted else 'none' }")
    print(f"  reads in accepted lengths: {acc_reads:,} ({pct:.1f}%)")
    if len(accepted) >= 3:
        off, n = Counter(accepted.values()).most_common(1)[0]
        print(f"  {n} of {len(accepted)} accepted lengths share offset {off}")

    fails, warns = [], []
    if not accepted:
        fails.append("no read length accepted as periodic")
    else:
        if pct < args.min_pct:
            fails.append(f"accepted lengths hold {pct:.1f}% of reads (minimum {args.min_pct:g}%)")
        if acc_reads < args.min_reads:
            warns.append(f"{acc_reads:,} reads in accepted lengths (fewer than {args.min_reads:,})")
        bad = [(k, v) for k, v in sorted(accepted.items()) if not plausible(k, v)]
        if bad:
            warns.append(f"implausible length/offset pairs: {bad}")
    verdict = "FAIL" if fails else ("WARN" if warns else "PASS")
    print(f"  verdict: {verdict}")
    for m in fails:
        print(f"    FAIL: {m}")
    for m in warns:
        print(f"    WARN: {m}")
    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("sample\tverdict\ttested_reads\taccepted_reads\tpct_accepted\taccepted_lengths\toffsets\tnotes\n")
            fh.write(f"{args.label}\t{verdict}\t{total}\t{acc_reads}\t{pct:.1f}\t"
                     f"{','.join(str(k) for k in sorted(accepted)) or '-'}\t"
                     f"{','.join(f'{k}:{accepted[k]}' for k in sorted(accepted)) or '-'}\t"
                     f"{'; '.join(fails + warns) or '-'}\n")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
