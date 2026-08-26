#!/usr/bin/env python3

"""Refold every sequence of an alignment under the consensus structure.

A reimplementation of ViennaRNA's ``refold.pl`` that needs nothing but the
ViennaRNA Python bindings: no Perl, and no ``RNAfold`` subprocess either,
because the constrained fold runs in-process through ``RNA.fold_compound``.

The constraint for a sequence is the consensus structure with

* every gap column removed,
* every pair to a gap opened,
* every pair the sequence cannot actually form opened, and
* every pair closing a loop shorter than ``--turn`` opened,

which is what ``refold.pl`` does and why its output is fed to ``RNAfold -C``.
Reading a dot plot (``RNAalifold -p``) keeps only pairs above ``--threshold``.

Deviations from ``refold.pl``, all of them cases where the Perl script fails
rather than produces a different answer:

* Sequence rows are read case-insensitively. ``refold.pl`` matches ``[A-Z-]+``
  and silently reads *nothing* from a lower-case alignment, then blocks on
  stdin; the workflow works around that by upper-casing in ``strip_aln``.
* A pair whose opening partner sits in column 0 is opened when its closing
  base is a gap. ``refold.pl`` guards that branch with ``$pt[$p] > 0``, leaves
  the ``(`` unmatched, and dies in its own pair-table check.
* ``--turn`` takes a value. ``refold.pl`` declares it as a boolean flag, so
  ``--turn 5`` sets the loop length to 1 and leaves ``5`` as a file argument.
* ``.``, ``~`` and ``_`` count as gaps alongside ``-``.

For any input ``refold.pl`` handles, the constraints are identical.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.tools.alignment_io import GAP_CHARACTERS, Alignment, parse_clustal_blocks

DEFAULT_THRESHOLD = 0.9
DEFAULT_TURN = 3

# The pairs refold.pl accepts: the six canonical ones plus their DNA spellings.
PAIRABLE = frozenset({"AU", "UA", "GC", "CG", "GU", "UG", "AT", "TA", "GT", "TG"})

UNPAIRED = -1


class RefoldError(RuntimeError):
    pass


@dataclass(frozen=True)
class RefoldRecord:
    """One sequence of the alignment, degapped, with its own constraint."""

    name: str
    sequence: str
    constraint: str


def pair_table(structure: str) -> list[int]:
    """Partner index per position, ``-1`` where unpaired. Indices start at 0."""
    table = [UNPAIRED] * len(structure)
    stack: list[int] = []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                raise RefoldError(f"Unbalanced brackets in consensus structure: {structure}")
            opened = stack.pop()
            table[index] = opened
            table[opened] = index
    if stack:
        raise RefoldError(f"Unbalanced brackets in consensus structure: {structure}")
    return table


def consensus_from_alifold(text: str) -> str:
    """The consensus structure from ``RNAalifold``'s own output."""
    for line in text.splitlines():
        stripped = line.strip()
        if stripped and set(stripped.split()[0]) <= set("(.)"):
            return stripped.split()[0]
    raise RefoldError("No consensus structure line found in the RNAalifold output")


def consensus_from_dotplot(text: str, length: int, threshold: float = DEFAULT_THRESHOLD) -> str:
    """The consensus structure from an ``RNAalifold -p`` dot plot.

    Only ``lbox`` entries at or above ``threshold`` become pairs. Fields are
    read from the end of the line, so both the coloured alifold dot plot
    (``0.00 1.00 hsb i j value lbox``) and a plain ``i j value lbox`` parse.
    """
    consensus = ["."] * length
    for line in text.splitlines():
        if not line.rstrip().endswith("lbox"):
            continue
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            i, j, value = int(fields[-4]), int(fields[-3]), float(fields[-2])
        except ValueError:
            continue
        if value < threshold:
            continue
        if not (1 <= i <= length and 1 <= j <= length):
            raise RefoldError(f"Dot plot pair {i}-{j} lies outside the {length}-column alignment")
        consensus[i - 1] = "("
        consensus[j - 1] = ")"
    return "".join(consensus)


def read_consensus(path: str | Path, length: int, threshold: float = DEFAULT_THRESHOLD) -> str:
    """Read the consensus from a dot plot or from RNAalifold's output."""
    text = Path(path).read_text(encoding="utf-8")
    if text.startswith("%!PS"):
        return consensus_from_dotplot(text, length, threshold)
    return consensus_from_alifold(text)


def constrain(sequence: str, consensus: str, turn: int = DEFAULT_TURN) -> tuple[str, str]:
    """Degap one aligned sequence and open every pair it cannot form."""
    if len(sequence) != len(consensus):
        raise RefoldError(
            f"Sequence of {len(sequence)} columns does not match a consensus of {len(consensus)}"
        )

    row = sequence.upper()
    constraint = list(consensus)
    table = pair_table(consensus)
    removed = "\x00"  # marks a gap column, dropped once every pair is resolved

    for position, char in enumerate(row):
        partner = table[position]
        if char in GAP_CHARACTERS:
            constraint[position] = removed
            if partner != UNPAIRED and constraint[partner] != removed:
                constraint[partner] = "."
        elif partner > position and row[position] + row[partner] not in PAIRABLE:
            constraint[position] = constraint[partner] = "."

    degapped = "".join(char for char in row if char not in GAP_CHARACTERS)
    kept = "".join(char for char in constraint if char != removed)

    # Loops shorter than `turn` cannot close, so those pairs would make the
    # constraint unsatisfiable rather than merely tight.
    kept_list = list(kept)
    table = pair_table(kept)
    for position, partner in enumerate(table):
        if position > partner:
            continue
        if partner - position - 1 < turn:
            kept_list[position] = kept_list[partner] = "."
    return degapped, "".join(kept_list)


def refold_alignment(alignment: Alignment, consensus: str, turn: int = DEFAULT_TURN) -> list[RefoldRecord]:
    records = []
    for name in alignment.order:
        sequence, constraint = constrain(alignment.seqs[name], consensus, turn)
        records.append(RefoldRecord(name=name, sequence=sequence, constraint=constraint))
    return records


def read_alignment(path: str | Path) -> Alignment:
    blocks = parse_clustal_blocks(Path(path).read_text(encoding="utf-8"))
    alignments = [block for block in blocks if block.order]
    if not alignments:
        raise RefoldError(f"No Clustal alignment found in {path}")
    if len(alignments) > 1:
        raise RefoldError(f"{path} holds {len(alignments)} alignments; refolding needs exactly one")
    return alignments[0]


def fold_constrained(sequence: str, constraint: str) -> tuple[str, float]:
    """MFE structure and energy under the constraint, as ``RNAfold -C`` gives them."""
    try:
        import RNA
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise RefoldError(
            "The ViennaRNA Python bindings are required to fold. Install ViennaRNA "
            "with its Python module (conda: viennarna, pip: ViennaRNA)."
        ) from error

    fold_compound = RNA.fold_compound(sequence)
    fold_compound.hc_add_from_db(constraint, RNA.CONSTRAINT_DB_DEFAULT)
    structure, energy = fold_compound.mfe()
    return structure, energy


def format_constraints(records: list[RefoldRecord]) -> str:
    """The stream ``refold.pl`` writes, for piping into ``RNAfold -C``."""
    lines = []
    for record in records:
        lines += [f"> {record.name}", record.sequence, record.constraint]
    return "\n".join(lines) + "\n" if lines else ""


def format_folded(records: list[RefoldRecord]) -> str:
    """What ``refold.pl | RNAfold --noPS -C`` writes, folded in-process."""
    lines = []
    for record in records:
        structure, energy = fold_constrained(record.sequence, record.constraint)
        lines += [f"> {record.name}", record.sequence, f"{structure} ({energy:6.2f})"]
    return "\n".join(lines) + "\n" if lines else ""


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Refold each sequence of an alignment under the consensus structure. "
            "Replaces refold.pl and the RNAfold -C it feeds."
        )
    )
    parser.add_argument("--alignment", required=True, help="Clustal alignment")
    parser.add_argument(
        "--consensus",
        required=True,
        help="RNAalifold output, or the dot plot from RNAalifold -p",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=DEFAULT_THRESHOLD,
        help="Dot plot only: keep pairs at or above this probability (default: %(default)s)",
    )
    parser.add_argument(
        "--turn",
        type=int,
        default=DEFAULT_TURN,
        help="Open pairs closing a loop shorter than this (default: %(default)s)",
    )
    parser.add_argument(
        "--constraints-only",
        action="store_true",
        help="Write refold.pl's constraint stream instead of folding it.",
    )
    parser.add_argument("--output", help="Write here instead of stdout")
    args = parser.parse_args()

    alignment = read_alignment(args.alignment)
    consensus = read_consensus(args.consensus, alignment.length, args.threshold)
    records = refold_alignment(alignment, consensus, args.turn)
    if not records:
        raise RefoldError(f"{args.alignment} holds no sequences to refold")

    text = format_constraints(records) if args.constraints_only else format_folded(records)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RefoldError as error:
        print(f"refold: {error}", file=sys.stderr)
        raise SystemExit(1) from error
