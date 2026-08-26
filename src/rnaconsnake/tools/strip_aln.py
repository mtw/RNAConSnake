#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import replace

from rnaconsnake.tools.stockholm_utils import (
    StockholmRecord,
    iter_stockholm_lines,
    parse_stockholm_records,
    subset_record,
)

FORMAT_MAP = {
    "C": "clustalw",
    "S": "stockholm",
    "M": "maf",
    "F": "fasta",
}

VERSION = "0.2"


def _uppercase_sequences(record: StockholmRecord) -> StockholmRecord:
    r"""Upper-case residues, leaving names and annotation lines untouched.

    This is the pipeline's one canonical case-normalisation point. The external
    toolchain does not agree about case: ``alifoldz.pl`` matches alignment rows
    with ``[A-Z\-.]+`` and ``refold.pl`` with ``[A-Z\-]+``, both without ``/i``,
    so a lower-case alignment (MAFFT output, for instance) parses as *zero*
    sequences in either -- while ``RNAz`` and ``RNAalifold`` accept any case.
    Mixed behaviour like that produces tool-dependent silent corruption, so
    everything downstream of this step gets upper case.

    ``generated_files/orig/`` still holds the verbatim split alignment.
    """
    return replace(record, seqs={name: seq.upper() for name, seq in record.seqs.items()})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove redundant sequences from alignments.")
    parser.add_argument("-a", "--aln", required=True, help="Input alignment")
    parser.add_argument("-f", "--format", required=True, help="Alignment format code")
    parser.add_argument(
        "--nosingle", action="store_true", help="Keep all sequences if only one unique sequence remains"
    )
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.version:
        print(f"strip_aln {VERSION}")
        return 0

    fmt = FORMAT_MAP.get(args.format)
    if fmt != "stockholm":
        print("Only Stockholm format is currently supported by the Python replacement.", file=sys.stderr)
        return 2

    records = parse_stockholm_records(args.aln)
    if not records:
        return 0

    record = _uppercase_sequences(records[0])
    keep: list[str] = []
    seen: set[str] = set()
    for name in record.seq_order:
        seq = record.seqs[name]
        if seq in seen:
            continue
        keep.append(name)
        seen.add(seq)

    if len(keep) == 1 and args.nosingle:
        subset = subset_record(record, record.seq_order)
    else:
        subset = subset_record(record, keep)

    for line in iter_stockholm_lines(subset):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
