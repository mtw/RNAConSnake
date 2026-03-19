#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

VERSION = "0.1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Split multi-Stockholm alignments into single-record files.")
    parser.add_argument("-a", "--aln", required=True, help="Input multi-Stockholm alignment")
    parser.add_argument("--accession", action="store_true", help="Prefer #=GF AC over #=GF ID for filenames")
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    return parser.parse_args()


def record_name(lines: list[str], prefer_accession: bool) -> str:
    record_id = None
    accession = None
    for line in lines:
        if line.startswith("#=GF AC "):
            accession = line.split(None, 2)[2]
        if line.startswith("#=GF ID "):
            record_id = line.split(None, 2)[2]
    if prefer_accession and accession:
        return accession
    if accession:
        return accession
    if record_id:
        return record_id
    raise ValueError("Record is missing both #=GF ID and #=GF AC")


def main() -> int:
    args = parse_args()
    if args.version:
        print(f"split_stockholm {VERSION}")
        return 0

    infile = Path(args.aln)
    if not infile.is_file():
        print("Could not find input file provided via --aln|-a option", file=sys.stderr)
        return 2

    current: list[str] = []
    with infile.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "//":
                if current:
                    out_path = Path(f"{record_name(current, args.accession)}.stk")
                    with out_path.open("w", encoding="utf-8") as out_handle:
                        for item in current:
                            out_handle.write(f"{item}\n")
                        out_handle.write("//\n")
                    current = []
                continue
            current.append(line)

    if current:
        out_path = Path(f"{record_name(current, args.accession)}.stk")
        with out_path.open("w", encoding="utf-8") as out_handle:
            for item in current:
                out_handle.write(f"{item}\n")
            out_handle.write("//\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
