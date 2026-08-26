#!/usr/bin/env python3

"""Report which regions of an input alignment can actually be screened.

A local RNA-structure screen is blind in two situations, and neither announces
itself: where most sequences are gaps, ``RNAalifold`` has no consensus to score
and reports nothing; and where mean pairwise identity is very low, ``RNAz``
loses discriminating power. Both look identical to a genuine absence of
structure -- the screen simply returns fewer candidates.

This report makes that explicit, so "no candidates here" can be told apart from
"this region was never screenable".
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from rnaconsnake.tools.alignment_io import (
    GAP_CHARACTERS,
    Alignment,
    parse_clustal_blocks,
    read_stockholm_alignment,
)
from rnaconsnake.workflow_helpers import write_json

# Below this fraction of sequences carrying a residue, a column contributes
# little but noise to the consensus fold.
DEFAULT_MIN_OCCUPANCY = 0.5
# RNAz's SVM is trained on alignments well above this; below it, class
# probabilities are not trustworthy.
DEFAULT_MIN_MPI = 0.55
DEFAULT_WINDOW = 40


@dataclass(frozen=True)
class RegionReport:
    start: int
    end: int
    occupancy: float
    mean_pairwise_identity: float
    screenable: bool
    reduced_power: bool
    reasons: tuple[str, ...]


def column_occupancy(alignment: Alignment, index: int) -> float:
    if not alignment.n_seq:
        return 0.0
    residues = sum(1 for name in alignment.order if alignment.seqs[name][index] not in GAP_CHARACTERS)
    return residues / alignment.n_seq


def window_identity(alignment: Alignment, start: int, end: int) -> float:
    if alignment.n_seq < 2:
        return 1.0
    matches = totals = 0
    for left, right in combinations(alignment.order, 2):
        a = alignment.seqs[left][start:end]
        b = alignment.seqs[right][start:end]
        for x, y in zip(a, b, strict=True):
            if x in GAP_CHARACTERS and y in GAP_CHARACTERS:
                continue
            totals += 1
            if x.upper() == y.upper():
                matches += 1
    return matches / totals if totals else 0.0


def profile(
    alignment: Alignment,
    window: int = DEFAULT_WINDOW,
    min_occupancy: float = DEFAULT_MIN_OCCUPANCY,
    min_identity: float = DEFAULT_MIN_MPI,
) -> list[RegionReport]:
    occupancies = [column_occupancy(alignment, i) for i in range(alignment.length)]
    regions: list[RegionReport] = []
    for start in range(0, alignment.length, window):
        end = min(start + window, alignment.length)
        occupancy = sum(occupancies[start:end]) / (end - start)
        identity = window_identity(alignment, start, end)
        # Only occupancy is disqualifying: with most sequences gapped there is
        # no consensus to fold, so nothing can be reported. Low identity merely
        # costs RNAz power -- strong elements are still found below the
        # threshold -- so it is a caveat, not a verdict.
        reasons: list[str] = []
        screenable = occupancy >= min_occupancy
        if not screenable:
            reasons.append(f"occupancy {occupancy:.2f} < {min_occupancy}")
        reduced_power = identity < min_identity
        if reduced_power:
            reasons.append(f"mean pairwise identity {identity:.2f} < {min_identity}")
        regions.append(
            RegionReport(
                start=start + 1,
                end=end,
                occupancy=round(occupancy, 4),
                mean_pairwise_identity=round(identity, 4),
                screenable=screenable,
                reduced_power=reduced_power,
                reasons=tuple(reasons),
            )
        )
    return regions


def summarize(regions: list[RegionReport], alignment: Alignment) -> dict:
    blind = [region for region in regions if not region.screenable]
    weak = [region for region in regions if region.screenable and region.reduced_power]
    blind_columns = sum(region.end - region.start + 1 for region in blind)
    weak_columns = sum(region.end - region.start + 1 for region in weak)
    warnings: list[str] = []
    if blind_columns:
        fraction = blind_columns / alignment.length
        warnings.append(
            f"{blind_columns} of {alignment.length} alignment columns ({fraction:.0%}) are too "
            "gappy to fold a consensus, so no candidate can be reported there. Absence of "
            "candidates in these regions is not evidence of absence of structure."
        )
        for start, end in _contiguous_runs(blind)[:5]:
            warnings.append(f"unscreenable (gappy): columns {start}-{end}")
    if weak_columns:
        fraction = weak_columns / alignment.length
        warnings.append(
            f"{weak_columns} columns ({fraction:.0%}) are below the identity range where RNAz "
            "discriminates well; candidates there are still reported but with reduced power."
        )
    return {
        "n_seq": alignment.n_seq,
        "length": alignment.length,
        "screenable_columns": alignment.length - blind_columns,
        "blind_columns": blind_columns,
        "reduced_power_columns": weak_columns,
        "warnings": warnings,
        "regions": [
            {
                "start": region.start,
                "end": region.end,
                "occupancy": region.occupancy,
                "mean_pairwise_identity": region.mean_pairwise_identity,
                "screenable": region.screenable,
                "reduced_power": region.reduced_power,
                "reasons": list(region.reasons),
            }
            for region in regions
        ],
    }


def _contiguous_runs(regions: list[RegionReport]) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    for region in regions:
        if runs and region.start == runs[-1][1] + 1:
            runs[-1] = (runs[-1][0], region.end)
        else:
            runs.append((region.start, region.end))
    return runs


def render_text(payload: dict) -> str:
    lines = [
        "# RNAcs alignment screenability report",
        f"# sequences: {payload['n_seq']}   columns: {payload['length']}",
        f"# screenable columns: {payload['screenable_columns']}   "
        f"unscreenable: {payload['blind_columns']}   "
        f"reduced power: {payload['reduced_power_columns']}",
        "#",
    ]
    for warning in payload["warnings"]:
        lines.append(f"# WARNING\t{warning}")
    lines.append("start\tend\toccupancy\tmean_pairwise_identity\tscreenable\treduced_power\treasons")
    for region in payload["regions"]:
        lines.append(
            "\t".join(
                [
                    str(region["start"]),
                    str(region["end"]),
                    f"{region['occupancy']:.4f}",
                    f"{region['mean_pairwise_identity']:.4f}",
                    "yes" if region["screenable"] else "no",
                    "yes" if region["reduced_power"] else "no",
                    "; ".join(region["reasons"]),
                ]
            )
        )
    return "\n".join(lines) + "\n"


def load_alignment(path: str | Path) -> Alignment:
    if Path(path).suffix.lower() == ".aln":
        blocks = parse_clustal_blocks(Path(path).read_text(encoding="utf-8"))
        if len(blocks) != 1:
            raise ValueError(f"Expected one Clustal alignment in {path}, found {len(blocks)}")
        return blocks[0]
    return read_stockholm_alignment(path)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report which regions of an alignment can be screened for structure."
    )
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--output", required=True, help="TSV report.")
    parser.add_argument("--metadata", help="Optional JSON copy of the same report.")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    parser.add_argument("--min-occupancy", type=float, default=DEFAULT_MIN_OCCUPANCY)
    parser.add_argument("--min-identity", type=float, default=DEFAULT_MIN_MPI)
    args = parser.parse_args()

    alignment = load_alignment(args.alignment)
    payload = summarize(profile(alignment, args.window, args.min_occupancy, args.min_identity), alignment)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(render_text(payload), encoding="utf-8")
    if args.metadata:
        write_json(args.metadata, payload)
    for warning in payload["warnings"]:
        print(f"WARNING: {warning}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
