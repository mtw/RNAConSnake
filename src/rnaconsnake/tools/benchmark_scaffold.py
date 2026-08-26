#!/usr/bin/env python3

"""Draft a benchmark truth file from a reference structure and a screen.

The positive control needs curated coordinates for known elements, and those
coordinates are alignment specific: they only mean something relative to the
exact alignment the benchmark runs on. Deriving them by hand is tedious and
easy to get subtly wrong.

This scaffolds the mechanical half. It takes a structure-annotated alignment --
an ``#=GC SS_cons`` line, as produced by mLocARNA -- and decomposes it into
element-scale structural domains, then cross-references each against what a
completed screen actually reported. What it emits is a truth file with the
coordinates filled in and the *labels* left blank.

Attaching the labels is a judgement call about biology and stays with the
curator: this tool will not guess which domain is xrRNA1.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.tools.benchmark import TRUTH_COLUMNS
from rnaconsnake.tools.dereplicate import base_pairs
from rnaconsnake.tools.loci import parse_window_name
from rnaconsnake.tools.stockholm_utils import parse_stockholm_records


# Structural elements of the kind this screen targets are on the order of tens
# of nucleotides. Helices spanning much more than this are long-range pairings
# that bracket several elements rather than defining one.
DEFAULT_MAX_WIDTH = 150
# Lowered from 15 after a confirmed 14 nt element (the dengue-group sHP) was
# silently dropped. A spurious row costs a curator a glance; a missing one costs
# an element nobody knows to look for.
DEFAULT_MIN_WIDTH = 10
DEFAULT_MIN_PAIRS = 3

PLACEHOLDER_ID = "TBD"


@dataclass(frozen=True)
class Domain:
    start: int
    end: int
    n_pairs: int

    @property
    def width(self) -> int:
        return self.end - self.start + 1


def read_ss_cons(path: str | Path) -> tuple[str, int]:
    records = parse_stockholm_records(path)
    if not records:
        raise ValueError(f"No Stockholm records in {path}")
    record = records[0]
    for line in record.gc_lines:
        if line.startswith("#=GC SS_cons"):
            parts = line.split(None, 2)
            if len(parts) == 3:
                return parts[2], len(record.seq_order)
    raise ValueError(
        f"{path} has no '#=GC SS_cons' line. A reference structure is required; "
        "run the alignment through mLocARNA (or add the consensus by hand) first."
    )


def structural_domains(
    structure: str,
    max_width: int = DEFAULT_MAX_WIDTH,
    min_width: int = DEFAULT_MIN_WIDTH,
    min_pairs: int = DEFAULT_MIN_PAIRS,
) -> list[Domain]:
    """Maximal helices no wider than ``max_width``.

    Descends through long-range pairings -- which bracket several elements
    rather than defining one -- and emits the first enclosing helix that is on
    an element scale.
    """
    pairs = base_pairs(structure)
    partner = dict(pairs)

    def count_pairs(lo: int, hi: int) -> int:
        return sum(1 for a, b in pairs if lo <= a and b <= hi)

    def branches(lo: int, hi: int) -> list[tuple[int, int]]:
        found, index = [], lo
        while index <= hi:
            if index in partner and partner[index] <= hi:
                found.append((index, partner[index]))
                index = partner[index] + 1
            else:
                index += 1
        return found

    domains: list[Domain] = []

    def visit(lo: int, hi: int) -> None:
        width = hi - lo + 1
        if width <= max_width:
            n_pairs = count_pairs(lo, hi)
            if width >= min_width and n_pairs >= min_pairs:
                domains.append(Domain(start=lo + 1, end=hi + 1, n_pairs=n_pairs))
            return
        for a, b in branches(lo + 1, hi - 1):
            visit(a, b)

    for a, b in branches(0, len(structure) - 1):
        visit(a, b)
    return sorted(domains, key=lambda domain: (domain.start, domain.end))


def read_reported_loci(paths: list[str | Path]) -> list[dict[str, str]]:
    loci: list[dict[str, str]] = []
    for path in paths:
        with open(path, encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("wbn"):
                    loci.append(row)
    return loci


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def _float(value: str | None) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def match_loci(domain: Domain, loci: list[dict[str, str]]) -> list[tuple[int, dict[str, str]]]:
    matches = []
    for row in loci:
        try:
            start = int(row["locus_start"])
            end = int(row["locus_end"])
        except (KeyError, TypeError, ValueError):
            coords = parse_window_name(row["wbn"])
            start, end = coords.start, coords.end
        columns = _overlap(domain.start, domain.end, start, end)
        if columns:
            matches.append((columns, row))
    return sorted(matches, key=lambda item: -item[0])


def scaffold_rows(
    domains: list[Domain],
    loci: list[dict[str, str]],
    alignment_id: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    matched_loci: set[str] = set()

    for index, domain in enumerate(domains, start=1):
        matches = match_loci(domain, loci)
        best = matches[0][1] if matches else None
        # Every overlapping locus counts as accounted for, not just the best
        # one: the same element is reported once per window length, and those
        # siblings are not unexplained findings.
        for _, row in matches:
            matched_loci.add(row.get("locus_id") or row["wbn"])
        note = [f"consensus domain, {domain.width} nt, {domain.n_pairs} bp"]
        if best is None:
            note.append("NOT REPORTED by the screen")
        else:
            probability = _float(best.get("rnazprob"))
            note.append(
                f"best overlapping locus {best.get('locus_id') or best['wbn']} "
                f"({matches[0][0]} nt overlap, rnazprob="
                f"{probability if probability is not None else 'NA'}, "
                f"alifoldz={best.get('alifoldzscore', 'NA')})"
            )
        rows.append(
            {
                "element_id": f"{PLACEHOLDER_ID}_{index:02d}",
                "element_class": PLACEHOLDER_ID,
                "alignment": alignment_id,
                "start": str(domain.start),
                "end": str(domain.end),
                "notes": "; ".join(note),
            }
        )

    # Loci the reference structure does not mark: the screen saw something the
    # consensus does not describe. Worth a decision, so they are listed too.
    for row in loci:
        identifier = row.get("locus_id") or row["wbn"]
        if identifier in matched_loci:
            continue
        probability = _float(row.get("rnazprob"))
        if probability is None or probability < 0.9:
            continue
        try:
            start = int(row["locus_start"])
            end = int(row["locus_end"])
        except (KeyError, TypeError, ValueError):
            coords = parse_window_name(row["wbn"])
            start, end = coords.start, coords.end
        rows.append(
            {
                "element_id": f"{PLACEHOLDER_ID}_extra_{identifier}",
                "element_class": PLACEHOLDER_ID,
                "alignment": alignment_id,
                "start": str(start),
                "end": str(end),
                "notes": (
                    f"reported by the screen (rnazprob={probability}) but not marked in the "
                    "reference structure; keep only if it is a real element"
                ),
            }
        )
    return rows


def write_scaffold(rows: list[dict[str, str]], path: str | Path, alignment_id: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(f"# Benchmark truth scaffold for alignment '{alignment_id}'.\n")
        handle.write("#\n")
        handle.write("# Coordinates are 1-based, inclusive ALIGNMENT COLUMNS and are filled in\n")
        handle.write("# from the reference consensus structure. The labels are NOT: replace every\n")
        handle.write(f"# '{PLACEHOLDER_ID}' element_id and element_class with the real element\n")
        handle.write("# name and class, drop any row that is not a real element, and merge or\n")
        handle.write("# split spans where the consensus decomposition disagrees with the biology.\n")
        handle.write("#\n")
        handle.write("# Rows still containing TBD are refused by the benchmark; see\n")
        handle.write("# resources/benchmark/README.md.\n")
        writer = csv.DictWriter(
            handle, fieldnames=TRUTH_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Draft a benchmark truth file from a reference structure and a screen."
    )
    parser.add_argument("--alignment", required=True, help="Stockholm alignment with #=GC SS_cons.")
    parser.add_argument("--alignment-id", required=True, help="Identifier for the truth file.")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--loci",
        action="append",
        default=[],
        help="RNAConSnake.nr.csv from a completed run. Repeatable.",
    )
    parser.add_argument("--max-width", type=int, default=DEFAULT_MAX_WIDTH)
    parser.add_argument("--min-width", type=int, default=DEFAULT_MIN_WIDTH)
    parser.add_argument("--min-pairs", type=int, default=DEFAULT_MIN_PAIRS)
    args = parser.parse_args()

    structure, n_seq = read_ss_cons(args.alignment)
    domains = structural_domains(structure, args.max_width, args.min_width, args.min_pairs)
    loci = read_reported_loci(args.loci)
    rows = scaffold_rows(domains, loci, args.alignment_id)
    write_scaffold(rows, args.output, args.alignment_id)

    # Lowering max_width does not improve which elements are captured -- a wide
    # domain already contains them -- but a curator cannot recognise a dumbbell
    # inside a 131 nt fold. Flag the wide ones for decomposition instead.
    oversized = [d for d in domains if d.width > 100]
    if oversized:
        print(
            f"note: {len(oversized)} domain(s) exceed 100 nt "
            f"({', '.join(f'{d.start}-{d.end} ({d.width}nt)' for d in oversized[:6])}). "
            "These may bracket several elements. Re-run with a smaller --max-width to "
            "decompose them, and fold the sub-spans with rnaconsnake.tools.fold_region."
        )
    small = [
        d for d in domains if d.width < 20
    ]
    if small:
        print(
            f"note: {len(small)} domain(s) are under 20 nt "
            f"({', '.join(f'{d.start}-{d.end}' for d in small[:6])}). Small hairpins are easy "
            "to miss and easy to invent; check each against the plots.",
        )
    missed = sum(1 for row in rows if "NOT REPORTED" in row["notes"])
    extra = sum(1 for row in rows if "not marked in the reference" in row["notes"])
    print(
        f"{n_seq} sequences, {len(structure)} columns -> {len(domains)} consensus domains "
        f"({len(domains) - missed} overlapped by a reported locus, {missed} missed), "
        f"plus {extra} reported loci absent from the reference -> {args.output}"
    )
    print("Labels are placeholders: edit element_id and element_class before use.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
