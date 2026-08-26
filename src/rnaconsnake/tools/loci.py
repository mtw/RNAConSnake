#!/usr/bin/env python3

"""Alignment-coordinate primitives shared by de-replication and calibration.

RNALalifold emits overlapping windows, so raw window counts are neither
independent nor free of redundancy. The clustering that turns windows into
loci lives in :mod:`rnaconsnake.tools.dereplicate`; this module holds the
coordinate handling both it and :mod:`rnaconsnake.tools.calibration` need.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# RNALalifold candidate basenames look like ``RC_150_0001_aln_2294_2375``:
# prefix, window length, block index, then the alignment coordinates.
WINDOW_NAME = re.compile(r"^(?P<prefix>.+)_aln_(?P<start>\d+)_(?P<end>\d+)$")


@dataclass(frozen=True)
class WindowCoords:
    prefix: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_window_name(name: str) -> WindowCoords:
    match = WINDOW_NAME.match(name)
    if not match:
        raise ValueError(f"Cannot parse alignment coordinates from candidate name {name!r}")
    start = int(match.group("start"))
    end = int(match.group("end"))
    if end < start:
        start, end = end, start
    return WindowCoords(prefix=match.group("prefix"), start=start, end=end)


def overlap(first: WindowCoords, second: WindowCoords) -> int:
    return min(first.end, second.end) - max(first.start, second.start) + 1


def collapse_ratio(n_windows: int, n_loci: int) -> float:
    """Windows per locus; ``0.0`` when there is nothing to collapse."""
    if n_loci <= 0:
        return 0.0
    return n_windows / n_loci
