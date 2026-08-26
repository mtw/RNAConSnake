#!/usr/bin/env python3

"""Positive control: recovery of known elements from a curated truth file.

The calibrated screen is only interpretable if it demonstrably *finds* known
structure.  This module compares the calibrated real-arm loci
(``results/calibration/qvalues.tsv``) against a hand-curated truth file of
known element coordinates and reports, per element, whether it was recovered
and at what q-value.

The truth file is repository data, not a generated artefact: see
``resources/benchmark/``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

TRUTH_COLUMNS = ["element_id", "element_class", "alignment", "start", "end", "notes"]

RECOVERY_COLUMNS = [
    "reciprocal_overlap",
    "element_id",
    "element_class",
    "alignment",
    "element_start",
    "element_end",
    "recovered",
    "best_locus",
    "locus_wlen",
    "locus_start",
    "locus_end",
    "overlap_columns",
    "overlap_fraction",
    "rnazprob",
    "alifoldzscore",
    "rscape_covary_count",
    "q_rnaz",
    "q_alifoldz",
    "q_cascade",
    "cascade_pass",
    "notes",
]

UNCURATED = {"", "TBD", "NA", "?"}

# benchmark_scaffold emits coordinate-complete rows whose labels are still
# placeholders. Those coordinates make `curated` true, so the labels are what
# must be checked -- otherwise an unedited scaffold runs as if it were curated
# and produces meaningless recovery numbers.
PLACEHOLDER_PREFIX = "TBD"

NA = "NA"


class BenchmarkError(RuntimeError):
    pass


@dataclass(frozen=True)
class TruthElement:
    element_id: str
    element_class: str
    alignment: str
    start: int | None
    end: int | None
    notes: str

    @property
    def curated(self) -> bool:
        return self.start is not None and self.end is not None

    @property
    def length(self) -> int:
        if not self.curated:
            return 0
        assert self.start is not None and self.end is not None
        return self.end - self.start + 1


def _coordinate(value: str | None) -> int | None:
    text = (value or "").strip()
    if text.upper() in UNCURATED:
        return None
    try:
        return int(text)
    except ValueError as error:
        raise BenchmarkError(f"Non-integer benchmark coordinate {value!r}") from error


def read_truth(path: str | Path, alignment: str | None = None) -> list[TruthElement]:
    rows: list[TruthElement] = []
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    missing = [column for column in TRUTH_COLUMNS if column not in (reader.fieldnames or [])]
    if missing:
        raise BenchmarkError(f"Truth file {path} is missing required columns: {', '.join(missing)}")
    for row in reader:
        element_alignment = (row.get("alignment") or "").strip()
        if alignment and element_alignment and element_alignment != alignment:
            continue
        rows.append(
            TruthElement(
                element_id=(row.get("element_id") or "").strip(),
                element_class=(row.get("element_class") or "").strip(),
                alignment=element_alignment,
                start=_coordinate(row.get("start")),
                end=_coordinate(row.get("end")),
                notes=(row.get("notes") or "").strip(),
            )
        )
    if not rows:
        raise BenchmarkError(f"Truth file {path} contains no elements for alignment {alignment!r}")
    return rows


def _read_locus_table(path: str | Path) -> list[dict[str, str]]:
    """A null arm's loci, from RNAConSnake.nr.csv."""
    with open(path, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [
        {
            "locus_id": row.get("locus_id", row.get("wbn", "")),
            "start": row.get("locus_start", row.get("start", "")),
            "end": row.get("locus_end", row.get("end", "")),
            "rnazprob": row.get("rnazprob", ""),
            "alifoldzscore": row.get("alifoldzscore", ""),
        }
        for row in rows
    ]


def read_qvalues(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def evaluate(
    truth: list[TruthElement],
    loci: list[dict[str, str]],
    min_overlap_fraction: float,
    allow_uncurated: bool,
) -> list[dict[str, object]]:
    placeholders = [
        element.element_id
        for element in truth
        if element.element_id.upper().startswith(PLACEHOLDER_PREFIX)
        or element.element_class.upper().startswith(PLACEHOLDER_PREFIX)
    ]
    if placeholders and not allow_uncurated:
        raise BenchmarkError(
            "Benchmark truth file still has placeholder labels for: "
            + ", ".join(placeholders[:10])
            + (" ..." if len(placeholders) > 10 else "")
            + ". This looks like an unedited benchmark_scaffold output: replace every "
            "TBD element_id and element_class with the real element name and class."
        )

    uncurated = [element.element_id for element in truth if not element.curated]
    if uncurated and not allow_uncurated:
        raise BenchmarkError(
            "Benchmark truth file still has uncurated coordinates for: "
            + ", ".join(uncurated)
            + ". Fill in resources/benchmark/ (see its README) or pass --allow-uncurated "
            "to emit them as 'uncurated' rows."
        )

    results: list[dict[str, object]] = []
    for element in truth:
        row: dict[str, object] = dict.fromkeys(RECOVERY_COLUMNS, NA)
        row.update(
            {
                "element_id": element.element_id,
                "element_class": element.element_class,
                "alignment": element.alignment or NA,
                "element_start": element.start if element.curated else NA,
                "element_end": element.end if element.curated else NA,
                "notes": element.notes,
            }
        )
        if not element.curated:
            row["recovered"] = "uncurated"
            results.append(row)
            continue

        assert element.start is not None and element.end is not None
        best: tuple[int, dict[str, str]] | None = None
        for locus in loci:
            try:
                locus_start = int(locus["start"])
                locus_end = int(locus["end"])
            except (KeyError, TypeError, ValueError):
                continue
            columns = _overlap(element.start, element.end, locus_start, locus_end)
            if columns <= 0:
                continue
            if best is None or columns > best[0]:
                best = (columns, locus)

        if best is None:
            row["recovered"] = "no"
            row["overlap_columns"] = 0
            row["overlap_fraction"] = "0"
            results.append(row)
            continue

        columns, locus = best
        fraction = columns / element.length if element.length else 0.0
        try:
            locus_length = int(locus["end"]) - int(locus["start"]) + 1
        except (KeyError, TypeError, ValueError):
            locus_length = 0
        reciprocal = min(fraction, columns / locus_length) if locus_length else 0.0
        row.update(
            {
                "recovered": "yes" if fraction >= min_overlap_fraction else "no",
                "best_locus": locus.get("locus_id", NA),
                "locus_wlen": locus.get("wlen", NA),
                "locus_start": locus.get("start", NA),
                "locus_end": locus.get("end", NA),
                "overlap_columns": columns,
                "overlap_fraction": f"{fraction:.4f}",
                # How much of the *locus* the element accounts for. A locus far
                # larger than the element still scores 1.0 on overlap_fraction,
                # which is how a screen reporting one huge locus can look
                # perfectly sensitive.
                "reciprocal_overlap": f"{reciprocal:.4f}",
                "rnazprob": locus.get("rnazprob", NA),
                "alifoldzscore": locus.get("alifoldzscore", NA),
                "rscape_covary_count": locus.get("rscape_covary_count", NA),
                "q_rnaz": locus.get("q_rnaz", NA),
                "q_alifoldz": locus.get("q_alifoldz", NA),
                "q_cascade": locus.get("q_cascade", NA),
                "cascade_pass": locus.get("cascade_pass", NA),
            }
        )
        results.append(row)
    return results


def null_baseline(
    truth: list[TruthElement],
    null_loci: list[list[dict[str, str]]],
    min_overlap_fraction: float,
) -> dict[str, float]:
    """How many elements a *null* arm's loci would also 'recover'.

    Reported candidates typically cover much of an alignment, so an overlap
    test can be satisfied by chance. Without this baseline a recovery count
    says nothing about whether the screen found real structure.
    """
    counts = [
        sum(
            1
            for row in evaluate(truth, loci, min_overlap_fraction, allow_uncurated=True)
            if row.get("recovered") == "yes"
        )
        for loci in null_loci
    ]
    if not counts:
        return {}
    return {
        "arms": len(counts),
        "mean": sum(counts) / len(counts),
        "min": min(counts),
        "max": max(counts),
    }


def write_recovery(
    results: list[dict[str, object]],
    path: str | Path,
    baseline: dict[str, float] | None = None,
) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        recovered = sum(1 for row in results if row.get("recovered") == "yes")
        curated = sum(1 for row in results if row.get("recovered") != "uncurated")
        handle.write(f"# elements\t{len(results)}\n")
        handle.write(f"# curated\t{curated}\n")
        handle.write(f"# recovered\t{recovered}\n")
        if baseline:
            handle.write(
                f"# null_baseline_recovered\tmean {baseline['mean']:.1f} "
                f"(min {baseline['min']}, max {baseline['max']}) over {baseline['arms']} null arms\n"
            )
            margin = recovered - baseline["mean"]
            handle.write(f"# recovery_margin_over_null\t{margin:+.1f}\n")
            if margin <= 0:
                handle.write(
                    "# WARNING\tthe null arms recover as many elements as the real arm; "
                    "this recovery count is not evidence of detection. Read the q-values instead.\n"
                )
            elif baseline["mean"] >= 0.8 * curated:
                handle.write(
                    "# WARNING\tnull arms recover most elements too, so the overlap test is "
                    "close to vacuous at this min-overlap-fraction. Raise it, or rely on "
                    "the q-values.\n"
                )
        writer = csv.DictWriter(handle, fieldnames=RECOVERY_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Score recovery of known elements against a curated truth file."
    )
    parser.add_argument("--truth", required=True)
    parser.add_argument("--qvalues", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--alignment", help="Restrict the truth file to one benchmark alignment.")
    parser.add_argument("--min-overlap-fraction", type=float, default=0.5)
    parser.add_argument(
        "--null-loci",
        action="append",
        nargs="+",
        default=[],
        metavar="TABLE",
        help=(
            "One null arm's locus tables (RNAConSnake.nr.csv), one per window length. "
            "Repeat the option once per arm. Adds a baseline showing how many elements "
            "a null arm also 'recovers', without which a recovery count is "
            "uninterpretable. Tables listed under a single --null-loci are pooled into "
            "one arm; splitting an arm's window lengths across several --null-loci "
            "would count each as its own arm and understate the baseline."
        ),
    )
    parser.add_argument(
        "--allow-uncurated",
        action="store_true",
        help="Emit rows for elements whose coordinates have not been curated yet.",
    )
    args = parser.parse_args()

    truth = read_truth(args.truth, args.alignment)
    loci = read_qvalues(args.qvalues)
    results = evaluate(truth, loci, args.min_overlap_fraction, args.allow_uncurated)

    null_loci = [[row for path in group for row in _read_locus_table(path)] for group in args.null_loci]
    baseline = null_baseline(truth, null_loci, args.min_overlap_fraction) if null_loci else None
    write_recovery(results, args.output, baseline)

    recovered = sum(1 for row in results if row.get("recovered") == "yes")
    curated = sum(1 for row in results if row.get("recovered") != "uncurated")
    print(f"Recovered {recovered}/{curated} curated benchmark elements -> {args.output}")
    if baseline:
        print(
            f"Null-arm baseline: {baseline['mean']:.1f}/{curated} over {baseline['arms']} arms "
            f"(margin {recovered - baseline['mean']:+.1f})"
        )
    if curated and recovered == 0:
        print("WARNING: no curated benchmark element was recovered", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
