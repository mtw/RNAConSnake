#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys

from rnaconsnake.tools.stockholm_utils import iter_stockholm_lines, parse_stockholm_records, subset_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove gap-only or high-gap sequences from Stockholm alignments."
    )
    parser.add_argument("-a", "--aln", required=True, help="Input alignment")
    parser.add_argument("-i", "--infmt", default="clustalw", help="Input alignment format")
    parser.add_argument("-o", "--outfmt", help="Output alignment format")
    parser.add_argument("-r", "--gapratio", type=float, default=1.0, help="Maximum allowed gap fraction")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    infmt = args.infmt.lower()
    outfmt = (args.outfmt or args.infmt).lower()
    if infmt != "stockholm" or outfmt != "stockholm":
        print("Only stockholm format is currently supported by the Python replacement.", file=sys.stderr)
        return 2
    if not (0 <= args.gapratio <= 1):
        print("gap ratio parameter must be 0 <= r <= 1", file=sys.stderr)
        return 2

    records = parse_stockholm_records(args.aln)
    if not records:
        return 0

    record = records[0]
    length = len(next(iter(record.seqs.values()))) if record.seqs else 0
    gap_string = "-" * length
    keep: list[str] = []

    for name in record.seq_order:
        seq = record.seqs[name]
        if seq == gap_string:
            continue
        gap_count = seq.count("-")
        if gap_count <= args.gapratio * length:
            keep.append(name)

    subset = subset_record(record, keep)
    for line in iter_stockholm_lines(subset):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
