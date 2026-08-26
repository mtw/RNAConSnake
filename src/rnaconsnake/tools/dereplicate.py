#!/usr/bin/env python3

"""Collapse RNALalifold's overlapping windows into non-redundant candidates.

RNALalifold reports every locally stable structure it finds, so a single real
element is typically reported many times over: once at close to its true
extent, and again as a series of shorter windows covering its individual
stable helices. Those fragments are themselves thermodynamically stable and
conserved, so they score well too, and the ranked candidate list ends up
dominated by pieces of the same element.

De-replication groups windows that describe the same element and reports one
representative for each, keeping every member on record so nothing is hidden.

Methods
-------
``containment`` (default)
    Two windows describe the same element when one interval is nested inside
    the other. This deliberately does *not* group merely overlapping windows:
    two elements that abut, or share a few columns, stay separate. Contrast
    single-linkage overlap clustering, which chains a run of pairwise
    overlapping windows into one span far longer than any real element.
``substructure``
    Stricter and structure-aware: a window is a fragment only when its
    consensus base pairs are (largely) a subset of the other's. More
    conservative, because a shorter window often re-folds into slightly shifted
    helices for the same physical element and so fails the subset test.
``overlap``
    Single-linkage clustering on shared columns, the classic RNAz-style locus
    merge. Prone to chaining.
``none``
    No de-replication; every window is its own candidate.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

from rnaconsnake.tools.loci import WindowCoords, overlap, parse_window_name
from rnaconsnake.workflow_helpers import SUMMARY_FIELDS, write_json


METHODS = ("containment", "substructure", "overlap", "none")

# How a locus picks the member window it is represented by.
#
# The choice matters more than it looks: RNALalifold reports sub-windows of a
# multi-lobed element, and a rule that prefers the best-scoring window will
# happily represent a dumbbell by one lobe -- which then fails the cascade that
# the whole element would have passed.
#
# Selecting by a score and then thresholding on that same score is a selection
# bias, but the null arm undergoes identical selection, so the empirical FDR
# stays calibrated. It is the *reported* score that is optimistic, not the FDR.
REPRESENTATIVE_RULES = ("best_rnaz", "widest", "best_alifoldz", "widest_of_top_half")

OPENING_BRACKETS = "(<[{"
CLOSING_BRACKETS = ")>]}"
BRACKET_PARTNER = dict(zip(CLOSING_BRACKETS, OPENING_BRACKETS))

NR_EXTRA_COLUMNS = [
    "locus_id",
    "locus_start",
    "locus_end",
    "n_windows",
    "members",
]

NR_COLUMNS = NR_EXTRA_COLUMNS + SUMMARY_FIELDS


def base_pairs(structure: str, offset: int = 0) -> frozenset[tuple[int, int]]:
    """Base pairs of a dot-bracket string, in absolute alignment columns.

    Unbalanced brackets are tolerated: RNAalifold consensus structures are
    well-formed, but a truncated or malformed record must not abort a run.
    """
    stacks: dict[str, list[int]] = {bracket: [] for bracket in OPENING_BRACKETS}
    pairs: set[tuple[int, int]] = set()
    for index, char in enumerate(structure or ""):
        if char in OPENING_BRACKETS:
            stacks[char].append(index)
        elif char in CLOSING_BRACKETS:
            stack = stacks[BRACKET_PARTNER[char]]
            if stack:
                pairs.add((stack.pop() + offset, index + offset))
    return frozenset(pairs)


@dataclass(frozen=True)
class Candidate:
    name: str
    start: int
    end: int
    pairs: frozenset[tuple[int, int]]
    record: dict[str, str]
    rank: tuple[float, ...]

    @property
    def coords(self) -> WindowCoords:
        return WindowCoords(prefix=self.name, start=self.start, end=self.end)

    @property
    def width(self) -> int:
        return self.end - self.start + 1


@dataclass
class Cluster:
    locus_id: str
    members: list[Candidate] = field(default_factory=list)
    rule: str = "best_rnaz"

    @property
    def start(self) -> int:
        return min(member.start for member in self.members)

    @property
    def end(self) -> int:
        return max(member.end for member in self.members)

    @property
    def representative(self) -> Candidate:
        return select_representative(self.members, self.rule)


def _numeric(value: str | None, default: float = float("-inf")) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return default if math.isnan(parsed) else parsed


def rank_key(record: dict[str, str]) -> tuple[float, ...]:
    """Ranking used to pick a cluster representative; larger is better.

    RNAz class probability first, then AlifoldZ (negated, since more negative
    is more significant), then covariation support and finally window width, so
    that a tie between a fragment and its parent resolves towards the parent.
    """
    return (
        _numeric(record.get("rnazprob")),
        -_numeric(record.get("alifoldzscore"), default=float("inf")),
        _numeric(record.get("maxcovarval")),
        _numeric(record.get("alilen")),
    )


def select_representative(members: list[Candidate], rule: str = "best_rnaz") -> Candidate:
    """Pick the window that represents a locus. Names break ties, so the choice
    is deterministic and identical in every arm."""
    if rule not in REPRESENTATIVE_RULES:
        raise ValueError(
            f"Unknown representative rule {rule!r}. Expected one of: "
            + ", ".join(REPRESENTATIVE_RULES)
        )
    if rule == "widest":
        return max(members, key=lambda m: (m.width, m.rank, m.name))
    if rule == "best_alifoldz":
        # More negative is better; missing sorts last.
        return min(
            members,
            key=lambda m: (
                _numeric(m.record.get("alifoldzscore"), default=float("inf")),
                -m.width,
                m.name,
            ),
        )
    if rule == "widest_of_top_half":
        ranked = sorted(members, key=lambda m: (m.rank, m.name), reverse=True)
        shortlist = ranked[: max(1, len(ranked) // 2)]
        return max(shortlist, key=lambda m: (m.width, m.rank, m.name))
    return max(members, key=lambda m: (m.rank, m.name))


def interval_contains(outer: Candidate, inner: Candidate) -> bool:
    return outer.start <= inner.start and inner.end <= outer.end


def pair_containment(inner: Candidate, outer: Candidate) -> float:
    """Fraction of ``inner``'s base pairs that ``outer`` also forms."""
    if not inner.pairs:
        return 0.0
    return len(inner.pairs & outer.pairs) / len(inner.pairs)


class _UnionFind:
    def __init__(self, size: int) -> None:
        self._parent = list(range(size))

    def find(self, item: int) -> int:
        while self._parent[item] != item:
            self._parent[item] = self._parent[self._parent[item]]
            item = self._parent[item]
        return item

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self._parent[right_root] = left_root


def _may_absorb(
    outer: Candidate, inner: Candidate, max_container_width: int, min_coverage: float
) -> bool:
    """May ``outer`` absorb the nested ``inner`` as a fragment of itself?

    RNALalifold's long windows can span several genuinely distinct elements: a
    200 nt window over three 60 nt hairpins is a scan artifact, not a parent,
    and letting it absorb them reports one locus where there are three.

    A window narrower than ``max_container_width`` is on the scale of a real
    element and may absorb anything nested in it. A wider one may only absorb
    windows that cover most of it -- so near-duplicate long windows still
    collapse into each other, while the small elements inside them do not.
    ``max_container_width`` of 0 disables the guard entirely.
    """
    if max_container_width <= 0 or outer.width <= max_container_width:
        return True
    return inner.width / outer.width >= min_coverage


def _linked(
    first: Candidate,
    second: Candidate,
    method: str,
    pair_threshold: float,
    min_overlap: int,
    max_container_width: int = 0,
    min_coverage: float = 0.8,
) -> bool:
    if method == "containment":
        return (
            interval_contains(first, second)
            and _may_absorb(first, second, max_container_width, min_coverage)
        ) or (
            interval_contains(second, first)
            and _may_absorb(second, first, max_container_width, min_coverage)
        )
    if method == "substructure":
        return (
            pair_containment(second, first) >= pair_threshold
            or pair_containment(first, second) >= pair_threshold
        )
    if method == "overlap":
        return overlap(first.coords, second.coords) >= min_overlap
    raise ValueError(f"Unknown de-replication method {method!r}")


def cluster_candidates(
    candidates: list[Candidate],
    method: str = "containment",
    pair_threshold: float = 0.9,
    min_overlap: int = 1,
    label: str = "locus",
    max_container_width: int = 0,
    min_coverage: float = 0.8,
    representative_rule: str = "best_rnaz",
) -> list[Cluster]:
    if method not in METHODS:
        raise ValueError(
            f"Unknown de-replication method {method!r}. Expected one of: " + ", ".join(METHODS)
        )
    ordered = sorted(candidates, key=lambda candidate: (candidate.start, candidate.end, candidate.name))
    if method == "none":
        return [
            Cluster(locus_id=f"{label}_{index:04d}", members=[candidate], rule=representative_rule)
            for index, candidate in enumerate(ordered, start=1)
        ]

    union = _UnionFind(len(ordered))
    for i, first in enumerate(ordered):
        for j in range(i + 1, len(ordered)):
            second = ordered[j]
            # Sorted by start, so once a window begins past the current one's
            # end it cannot be nested in it -- but a later window may still
            # contain it, so only the overlap method can stop early.
            if method == "overlap" and second.start > first.end:
                break
            if _linked(
                first,
                second,
                method,
                pair_threshold,
                min_overlap,
                max_container_width,
                min_coverage,
            ):
                union.union(i, j)

    grouped: dict[int, list[Candidate]] = {}
    for index, candidate in enumerate(ordered):
        grouped.setdefault(union.find(index), []).append(candidate)

    clusters = [
        Cluster(locus_id="", members=members, rule=representative_rule)
        for members in grouped.values()
    ]
    clusters.sort(key=lambda cluster: (cluster.start, cluster.end))
    for index, cluster in enumerate(clusters, start=1):
        cluster.locus_id = f"{label}_{index:04d}"
    return clusters


def candidates_from_records(records: list[dict[str, str]]) -> list[Candidate]:
    candidates: list[Candidate] = []
    for record in records:
        name = (record.get("wbn") or "").strip()
        if not name:
            continue
        coords = parse_window_name(name)
        candidates.append(
            Candidate(
                name=name,
                start=coords.start,
                end=coords.end,
                pairs=base_pairs(record.get("alifold_consstruc", ""), coords.start),
                record=record,
                rank=rank_key(record),
            )
        )
    return candidates


def read_summary_csv(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_non_redundant_csv(clusters: list[Cluster], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(clusters, key=lambda cluster: (-cluster.representative.rank[0], cluster.start))
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NR_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for cluster in ordered:
            representative = cluster.representative
            row = {field: representative.record.get(field, "") for field in SUMMARY_FIELDS}
            row.update(
                {
                    "locus_id": cluster.locus_id,
                    "locus_start": cluster.start,
                    "locus_end": cluster.end,
                    "n_windows": len(cluster.members),
                    "members": ";".join(
                        sorted(member.name for member in cluster.members)
                    ),
                }
            )
            writer.writerow(row)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Collapse overlapping RNALalifold windows into non-redundant candidates."
    )
    parser.add_argument("--input", required=True, help="Per-window summary CSV.")
    parser.add_argument("--output", required=True, help="Non-redundant candidate CSV.")
    parser.add_argument("--method", default="containment", choices=list(METHODS))
    parser.add_argument("--pair-containment", type=float, default=0.9)
    parser.add_argument("--min-overlap", type=int, default=1)
    parser.add_argument(
        "--container-min-coverage",
        type=float,
        default=0.8,
        help=(
            "A window wider than --max-container-width may only absorb windows covering "
            "at least this fraction of it, so near-duplicate long windows still collapse."
        ),
    )
    parser.add_argument(
        "--max-container-width",
        type=int,
        default=0,
        help=(
            "A window wider than this may not absorb others; it forms its own locus. "
            "Stops one long RNALalifold window from swallowing several distinct "
            "elements. 0 disables the guard."
        ),
    )
    parser.add_argument("--label", default="locus")
    parser.add_argument(
        "--representative", default="best_rnaz", choices=list(REPRESENTATIVE_RULES),
        help="How each locus picks the window that represents it.",
    )
    parser.add_argument(
        "--metadata",
        help="Optional JSON sidecar recording the method and counts, for consumers.",
    )
    args = parser.parse_args()

    candidates = candidates_from_records(read_summary_csv(args.input))
    clusters = cluster_candidates(
        candidates,
        method=args.method,
        pair_threshold=args.pair_containment,
        min_overlap=args.min_overlap,
        label=args.label,
        max_container_width=args.max_container_width,
        min_coverage=args.container_min_coverage,
        representative_rule=args.representative,
    )
    write_non_redundant_csv(clusters, args.output)
    if args.metadata:
        write_json(
            args.metadata,
            {
                "method": args.method,
                "pair_containment": args.pair_containment,
                "min_overlap": args.min_overlap,
                "max_container_width": args.max_container_width,
                "container_min_coverage": args.container_min_coverage,
                "representative_rule": args.representative,
                "label": args.label,
                "n_windows": len(candidates),
                "n_loci": len(clusters),
            },
        )
    print(
        f"de-replicated {len(candidates)} windows into {len(clusters)} candidates "
        f"({args.method}) -> {args.output}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
