#!/usr/bin/env python3

"""Measure how far detection survives as an alignment gets harder.

A screen that returns nothing has two possible explanations: there is nothing
there, or the alignment could not have shown it. The null arm bounds false
positives but says nothing about this floor, so "no candidates" on a divergent
alignment is not by itself evidence of absence.

This measures the floor directly. It degrades an alignment whose elements are
known -- dropping sequences, and preferring subsets at a target identity -- and
records at what point the known elements stop being recovered. A screen that
still finds them at the divergence of the alignment under test can claim
evidence of absence there; one that does not, cannot.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import random
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.tools.alignment_io import (
    Alignment,
    mean_pairwise_identity,
    read_stockholm_alignment,
)


@dataclass(frozen=True)
class Subset:
    names: tuple[str, ...]
    identity: float

    @property
    def n_seq(self) -> int:
        return len(self.names)


def subset_identity(alignment: Alignment, names: tuple[str, ...]) -> float:
    return mean_pairwise_identity(Alignment(order=list(names), seqs={n: alignment.seqs[n] for n in names}))


def enumerate_subsets(
    alignment: Alignment,
    n_seq: int,
    max_subsets: int = 200,
    seed: int = 0,
) -> list[Subset]:
    """Subsets of a given size, with their mean pairwise identity.

    Exhaustive while that is cheap, otherwise a deterministic random sample --
    the point is to span the identity range, not to be complete.
    """
    if n_seq < 2 or n_seq > alignment.n_seq:
        return []
    names = list(alignment.order)
    total = 1
    for k in range(n_seq):
        total = total * (len(names) - k) // (k + 1)
    if total <= max_subsets:
        combos = list(itertools.combinations(names, n_seq))
    else:
        rng = random.Random(seed)
        seen: set[tuple[str, ...]] = set()
        while len(seen) < max_subsets:
            seen.add(tuple(sorted(rng.sample(names, n_seq))))
        combos = sorted(seen)
    return sorted(
        (Subset(names=tuple(c), identity=round(subset_identity(alignment, c), 4)) for c in combos),
        key=lambda s: s.identity,
    )


def pick_by_identity(subsets: list[Subset], targets: list[float]) -> list[Subset]:
    """The subset closest to each target identity, without repeats."""
    chosen: list[Subset] = []
    used: set[tuple[str, ...]] = set()
    for target in targets:
        remaining = [s for s in subsets if s.names not in used]
        if not remaining:
            break
        best = min(remaining, key=lambda s: abs(s.identity - target))
        chosen.append(best)
        used.add(best.names)
    return chosen


def write_subset_alignment(alignment: Alignment, subset: Subset, path: Path, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        handle.write(f"#=GF ID {label}\n")
        handle.writelines(f"{name} {alignment.seqs[name]}\n" for name in subset.names)
        handle.write("//\n")


def plan(
    alignment: Alignment,
    sizes: list[int],
    targets: list[float],
    max_subsets: int,
    seed: int,
) -> list[dict]:
    rows: list[dict] = []
    for n_seq in sizes:
        subsets = enumerate_subsets(alignment, n_seq, max_subsets, seed)
        if not subsets:
            continue
        for subset in pick_by_identity(subsets, targets):
            rows.append(
                {
                    "n_seq": n_seq,
                    "identity": subset.identity,
                    "names": ";".join(subset.names),
                    "label": f"n{n_seq}_id{round(subset.identity * 100):02d}",
                }
            )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plan a sensitivity-envelope experiment over alignment subsets."
    )
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--output-dir", required=True, help="Where subset alignments are written.")
    parser.add_argument("--manifest", required=True, help="TSV describing the subsets.")
    parser.add_argument("--sizes", default="", help="Comma-separated subset sizes. Default: 2..n.")
    parser.add_argument(
        "--identities",
        default="0.4,0.5,0.6,0.7,0.8,0.9",
        help="Target mean pairwise identities to sample toward.",
    )
    parser.add_argument("--max-subsets", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261101)
    args = parser.parse_args()

    alignment = read_stockholm_alignment(args.alignment)
    sizes = (
        [int(x) for x in args.sizes.split(",") if x.strip()]
        if args.sizes
        else list(range(2, alignment.n_seq + 1))
    )
    targets = [float(x) for x in args.identities.split(",") if x.strip()]
    rows = plan(alignment, sizes, targets, args.max_subsets, args.seed)

    outdir = Path(args.output_dir)
    for row in rows:
        subset = Subset(names=tuple(row["names"].split(";")), identity=row["identity"])
        write_subset_alignment(alignment, subset, outdir / f"{row['label']}.stk", row["label"])

    Path(args.manifest).parent.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, "w", encoding="utf-8", newline="") as handle:
        handle.write(f"# sensitivity envelope subsets of {args.alignment}\n")
        handle.write(f"# source: {alignment.n_seq} sequences, MPI {mean_pairwise_identity(alignment):.4f}\n")
        writer = csv.DictWriter(
            handle,
            fieldnames=["label", "n_seq", "identity", "names"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in ["label", "n_seq", "identity", "names"]})
    print(
        f"{len(rows)} subsets spanning "
        f"MPI {min(r['identity'] for r in rows):.2f}-{max(r['identity'] for r in rows):.2f} "
        f"-> {outdir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
