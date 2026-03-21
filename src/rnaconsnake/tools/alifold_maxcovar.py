#!/usr/bin/env python3

from __future__ import annotations

import sys
from typing import Iterable


def compute(lines: Iterable[str]) -> tuple[int, str]:
    maxcovar = 0
    gapali = 0
    cov = {i: 0 for i in range(7)}

    for raw_line in lines:
        line = raw_line.rstrip("\n")
        if not line:
            continue
        stripped = line.lstrip()
        if not stripped or not stripped[0].isdigit():
            continue

        data = line.split()
        if not data:
            continue

        last = data[-1]
        if last.endswith("+") or last.endswith("-"):
            continue

        has_gapali = last.startswith("--:") and last[3:].isdigit()
        if has_gapali:
            gapali += 1

        try:
            mismatch_count = int(data[2])
        except (IndexError, ValueError):
            continue

        if mismatch_count <= 2:
            covar = len(data) - 5
            if has_gapali:
                covar -= 1
            if covar > maxcovar:
                maxcovar = covar
            cov[covar] = cov.get(covar, 0) + 1

    what = f"gapali_{gapali}" if gapali > 0 else str(cov.get(maxcovar, 0))
    return maxcovar, what


def main() -> int:
    maxcovar, what = compute(sys.stdin)
    sys.stdout.write(f"{maxcovar}\t{what}")
    return maxcovar


if __name__ == "__main__":
    raise SystemExit(main())
