#!/usr/bin/env python3

"""Alignment format helpers shared by the null-model and calibration tools.

The pipeline itself uses ``esl-reformat`` for Stockholm/Clustal conversion, but
the null arm has to read *multi-block* Clustal produced by SISSIz and
``rnazRandomizeAln.pl`` (one block per simulated replicate), which is outside
what ``esl-reformat`` handles.  Keeping the parsing here also lets the null
generator compute composition diagnostics without another external tool.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from rnaconsnake.tools.stockholm_utils import StockholmRecord, parse_stockholm_records

GAP_CHARACTERS = frozenset("-.~_")
CONSERVATION_CHARACTERS = frozenset("*:. ")


@dataclass(frozen=True)
class Alignment:
    order: list[str]
    seqs: dict[str, str]

    @property
    def length(self) -> int:
        return len(self.seqs[self.order[0]]) if self.order else 0

    @property
    def n_seq(self) -> int:
        return len(self.order)

    def column(self, index: int) -> str:
        return "".join(self.seqs[name][index] for name in self.order)


def _is_conservation_line(line: str) -> bool:
    return bool(line) and set(line) <= CONSERVATION_CHARACTERS


def parse_clustal_blocks(text: str) -> list[Alignment]:
    """Split Clustal text into one :class:`Alignment` per ``CLUSTAL`` header."""
    blocks: list[Alignment] = []
    order: list[str] = []
    seqs: dict[str, str] = {}

    def flush() -> None:
        if order:
            blocks.append(Alignment(order=list(order), seqs=dict(seqs)))
        order.clear()
        seqs.clear()

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.lstrip().upper().startswith("CLUSTAL"):
            flush()
            continue
        if not line.strip():
            continue
        if line[0].isspace() or _is_conservation_line(line.strip()):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        name, chunk = parts[0], parts[1]
        if name not in seqs:
            order.append(name)
            seqs[name] = ""
        seqs[name] += chunk

    flush()
    return blocks


def write_clustal(alignment: Alignment, path: str | Path, width: int = 60) -> None:
    name_width = max((len(name) for name in alignment.order), default=1)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("CLUSTAL W (1.81) multiple sequence alignment\n\n\n")
        for start in range(0, alignment.length, width):
            for name in alignment.order:
                chunk = alignment.seqs[name][start : start + width]
                handle.write(f"{name.ljust(name_width)} {chunk}\n")
            handle.write("\n")


def alignment_from_stockholm_record(record: StockholmRecord) -> Alignment:
    return Alignment(order=list(record.seq_order), seqs=dict(record.seqs))


def read_stockholm_alignment(path: str | Path) -> Alignment:
    records = parse_stockholm_records(path)
    if not records:
        raise ValueError(f"No Stockholm records found in {path}")
    if len(records) > 1:
        raise ValueError(f"Expected a single alignment in {path}, found {len(records)} Stockholm records")
    return alignment_from_stockholm_record(records[0])


def write_stockholm_alignment(
    alignment: Alignment,
    path: str | Path,
    identifier: str | None = None,
    extra_gf: list[str] | None = None,
) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        if identifier:
            handle.write(f"#=GF ID {identifier}\n")
        for line in extra_gf or []:
            handle.write(f"{line}\n")
        for name in alignment.order:
            handle.write(f"{name} {alignment.seqs[name]}\n")
        handle.write("//\n")


def gap_mask(alignment: Alignment) -> list[str]:
    """Per-sequence gap masks, one ``'-'``/``'x'`` string per sequence."""
    return [
        "".join("-" if char in GAP_CHARACTERS else "x" for char in alignment.seqs[name])
        for name in alignment.order
    ]


def mean_pairwise_identity(alignment: Alignment) -> float:
    if alignment.n_seq < 2:
        return 1.0
    totals = 0
    matches = 0
    for left, right in combinations(alignment.order, 2):
        a, b = alignment.seqs[left], alignment.seqs[right]
        for char_a, char_b in zip(a, b, strict=True):
            if char_a in GAP_CHARACTERS and char_b in GAP_CHARACTERS:
                continue
            totals += 1
            if char_a.upper() == char_b.upper():
                matches += 1
    return matches / totals if totals else 1.0


def base_composition(alignment: Alignment) -> dict[str, float]:
    counts: Counter[str] = Counter()
    for name in alignment.order:
        for char in alignment.seqs[name].upper():
            if char in GAP_CHARACTERS:
                continue
            counts[("U" if char == "T" else char)] += 1
    total = sum(counts.values())
    if not total:
        return {}
    return {base: round(count / total, 6) for base, count in sorted(counts.items())}


def dinucleotide_composition(alignment: Alignment) -> dict[str, float]:
    counts: Counter[str] = Counter()
    for name in alignment.order:
        ungapped = "".join(
            ("U" if char == "T" else char)
            for char in alignment.seqs[name].upper()
            if char not in GAP_CHARACTERS
        )
        for index in range(len(ungapped) - 1):
            counts[ungapped[index : index + 2]] += 1
    total = sum(counts.values())
    if not total:
        return {}
    return {pair: round(count / total, 6) for pair, count in sorted(counts.items())}


def gap_fraction(alignment: Alignment) -> float:
    total = alignment.length * alignment.n_seq
    if not total:
        return 0.0
    gaps = sum(1 for name in alignment.order for char in alignment.seqs[name] if char in GAP_CHARACTERS)
    return round(gaps / total, 6)


def alignment_diagnostics(alignment: Alignment) -> dict[str, object]:
    return {
        "n_seq": alignment.n_seq,
        "length": alignment.length,
        "gap_fraction": gap_fraction(alignment),
        "mean_pairwise_identity": round(mean_pairwise_identity(alignment), 6),
        "base_composition": base_composition(alignment),
        "dinucleotide_composition": dinucleotide_composition(alignment),
    }


def uses_rna_alphabet(alignment: Alignment) -> bool:
    for name in alignment.order:
        upper = alignment.seqs[name].upper()
        if "U" in upper:
            return True
        if "T" in upper:
            return False
    return True
