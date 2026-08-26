#!/usr/bin/env python3

"""Sweep the cascade thresholds against the null arm.

The RNAz and AlifoldZ cutoffs decide what a screen reports, but they are
usually inherited rather than measured. Once a null arm exists the question is
answerable: at each candidate threshold pair, how many real loci survive, how
many survive by chance, and what is the resulting FDR?

This reuses the per-arm summary tables a calibrated run already produced, so a
sweep costs seconds rather than another screen.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.tools.calibration import (
    LocusRecord,
    Thresholds,
    build_loci,
    cascade_survivors,
    read_summary_rows,
    rscape_was_evaluated,
)
from rnaconsnake.workflow_helpers import REAL_ARM


SWEEP_COLUMNS = [
    "rnaz_prob",
    "alifoldz",
    "real_survivors",
    "null_survivors_mean",
    "null_survivors_sd",
    "fdr",
    "recovered_reference",
    "reference_total",
]


@dataclass(frozen=True)
class SweepPoint:
    rnaz_prob: float
    alifoldz: float
    real_survivors: int
    null_mean: float
    null_sd: float
    fdr: float | None
    recovered: int | None
    reference_total: int | None


def _overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def read_reference_spans(path: str | Path) -> list[tuple[int, int]]:
    """Reference element spans, for the recovery column. Labels are ignored."""
    with open(path, encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    spans: list[tuple[int, int]] = []
    for row in csv.DictReader(lines, delimiter="\t"):
        try:
            spans.append((int(row["start"]), int(row["end"])))
        except (KeyError, TypeError, ValueError):
            continue
    return spans


def count_recovered(
    survivors: list[LocusRecord], spans: list[tuple[int, int]], min_fraction: float
) -> int:
    recovered = 0
    for start, end in spans:
        length = end - start + 1
        best = max(
            (_overlap(start, end, locus.start, locus.end) for locus in survivors), default=0
        )
        if length and best / length >= min_fraction:
            recovered += 1
    return recovered


def sweep(
    arm_inputs: dict[str, dict[int, str]],
    rnaz_grid: list[float],
    alifoldz_grid: list[float],
    base: Thresholds,
    reference_spans: list[tuple[int, int]] | None = None,
    min_fraction: float = 0.5,
) -> list[SweepPoint]:
    # Cluster once: the loci do not depend on the score thresholds, only the
    # survivor counts do.
    loci_by_arm: dict[str, list[LocusRecord]] = {}
    for arm, per_wlen in arm_inputs.items():
        loci: list[LocusRecord] = []
        for wlen, path in sorted(per_wlen.items()):
            loci.extend(
                build_loci(
                    read_summary_rows(path),
                    arm,
                    wlen,
                    base.dereplicate_method,
                    base.pair_containment,
                    base.locus_min_overlap,
                    base.max_container_width,
                    base.container_min_coverage,
                )
            )
        loci_by_arm[arm] = loci

    null_arms = sorted(arm for arm in loci_by_arm if arm != REAL_ARM)
    include_rscape = rscape_was_evaluated(loci_by_arm)

    points: list[SweepPoint] = []
    for rnaz_prob in rnaz_grid:
        for alifoldz in alifoldz_grid:
            thresholds = Thresholds(
                rnaz_prob=rnaz_prob,
                alifoldz=alifoldz,
                rscape_min_pairs=base.rscape_min_pairs,
                stage1_rnaz_prob=min(base.stage1_rnaz_prob, rnaz_prob),
                locus_min_overlap=base.locus_min_overlap,
                collapse_ratio_tolerance=base.collapse_ratio_tolerance,
                dereplicate_method=base.dereplicate_method,
                pair_containment=base.pair_containment,
                max_container_width=base.max_container_width,
                container_min_coverage=base.container_min_coverage,
            )
            real = cascade_survivors(loci_by_arm[REAL_ARM], thresholds, include_rscape)
            null_counts = [
                len(cascade_survivors(loci_by_arm[arm], thresholds, include_rscape))
                for arm in null_arms
            ]
            null_mean = statistics.fmean(null_counts) if null_counts else 0.0
            null_sd = statistics.pstdev(null_counts) if len(null_counts) > 1 else 0.0
            fdr = min(1.0, null_mean / len(real)) if real else None
            recovered = (
                count_recovered(real, reference_spans, min_fraction)
                if reference_spans is not None
                else None
            )
            points.append(
                SweepPoint(
                    rnaz_prob=rnaz_prob,
                    alifoldz=alifoldz,
                    real_survivors=len(real),
                    null_mean=round(null_mean, 4),
                    null_sd=round(null_sd, 4),
                    fdr=None if fdr is None else round(fdr, 6),
                    recovered=recovered,
                    reference_total=len(reference_spans) if reference_spans is not None else None,
                )
            )
    return points


def write_sweep(points: list[SweepPoint], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("# cascade threshold sweep against the null arm\n")
        handle.write("# recovered_reference counts reference spans overlapped by a survivor;\n")
        handle.write("# it measures sensitivity only and says nothing about precision.\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SWEEP_COLUMNS)
        for point in points:
            writer.writerow(
                [
                    point.rnaz_prob,
                    point.alifoldz,
                    point.real_survivors,
                    point.null_mean,
                    point.null_sd,
                    "NA" if point.fdr is None else point.fdr,
                    "NA" if point.recovered is None else point.recovered,
                    "NA" if point.reference_total is None else point.reference_total,
                ]
            )


def _parse_grid(text: str) -> list[float]:
    return [float(part) for part in text.split(",") if part.strip()]


def _parse_arm_input(token: str) -> tuple[str, int, str]:
    arm, wlen, path = token.split(":", 2)
    return arm, int(wlen), path


def main() -> int:
    parser = argparse.ArgumentParser(description="Sweep cascade thresholds against the null arm.")
    parser.add_argument("--arm-input", action="append", required=True, metavar="ARM:WLEN:PATH")
    parser.add_argument("--output", required=True)
    parser.add_argument("--rnaz-grid", default="0.5,0.6,0.7,0.8,0.9,0.95,0.99")
    parser.add_argument("--alifoldz-grid", default="0,-1,-1.5,-2,-2.5,-3")
    parser.add_argument("--reference", help="Truth or scaffold TSV; adds a recovery column.")
    parser.add_argument("--min-overlap-fraction", type=float, default=0.5)
    parser.add_argument("--dereplicate-method", default="containment")
    parser.add_argument("--max-container-width", type=int, default=120)
    args = parser.parse_args()

    arm_inputs: dict[str, dict[int, str]] = {}
    for token in args.arm_input:
        arm, wlen, path = _parse_arm_input(token)
        arm_inputs.setdefault(arm, {})[wlen] = path

    base = Thresholds(
        rnaz_prob=0.9,
        alifoldz=-2.0,
        rscape_min_pairs=1,
        stage1_rnaz_prob=0.5,
        locus_min_overlap=1,
        collapse_ratio_tolerance=0.2,
        dereplicate_method=args.dereplicate_method,
        pair_containment=0.9,
        max_container_width=args.max_container_width,
        container_min_coverage=0.8,
    )
    spans = read_reference_spans(args.reference) if args.reference else None
    points = sweep(
        arm_inputs,
        _parse_grid(args.rnaz_grid),
        _parse_grid(args.alifoldz_grid),
        base,
        spans,
        args.min_overlap_fraction,
    )
    write_sweep(points, args.output)
    print(f"{len(points)} threshold combinations -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
