#!/usr/bin/env python3

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


@dataclass
class StockholmRecord:
    gf_lines: list[str]
    gc_lines: list[str]
    seq_order: list[str]
    seqs: dict[str, str]
    gs_lines: list[tuple[str, str]]
    gr_lines: list[tuple[str, str]]
    other_lines: list[str]


def _split_seq_line(line: str) -> tuple[str, str]:
    parts = line.rstrip("\n").split(None, 1)
    if len(parts) != 2:
        raise ValueError(f"Malformed sequence line: {line!r}")
    return parts[0], parts[1]


def parse_stockholm_records(path: str | Path) -> list[StockholmRecord]:
    records: list[StockholmRecord] = []
    current = _new_record()

    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line == "//":
                if (
                    current["seq_order"]
                    or current["gf_lines"]
                    or current["gc_lines"]
                    or current["other_lines"]
                ):
                    records.append(_record_from_parts(current))
                current = _new_record()
                continue

            if line.startswith("#=GF "):
                current["gf_lines"].append(line)
                continue
            if line.startswith("#=GC "):
                _append_gc(current, line)
                continue
            if line.startswith("#=GS "):
                parts = line.split(None, 2)
                if len(parts) != 3:
                    raise ValueError(f"Malformed #=GS line: {line!r}")
                current["gs_lines"].append((parts[1], line))
                continue
            if line.startswith("#=GR "):
                parts = line.split(None, 2)
                if len(parts) != 3:
                    raise ValueError(f"Malformed #=GR line: {line!r}")
                current["gr_lines"].append((parts[1], line))
                continue
            if line.startswith("#"):
                current["other_lines"].append(line)
                continue

            name, seq = _split_seq_line(line)
            if name not in current["seqs"]:
                current["seq_order"].append(name)
                current["seqs"][name] = seq
            else:
                # Interleaved Stockholm repeats each name once per block.
                current["seqs"][name] += seq

    if current["seq_order"] or current["gf_lines"] or current["gc_lines"] or current["other_lines"]:
        records.append(_record_from_parts(current))
    return records


def write_stockholm_record(record: StockholmRecord, path: str | Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(f"{line}\n" for line in iter_stockholm_lines(record))


def iter_stockholm_lines(record: StockholmRecord) -> Iterable[str]:
    yield "# STOCKHOLM 1.0"
    for line in record.other_lines:
        if line != "# STOCKHOLM 1.0":
            yield line
    for line in record.gf_lines:
        yield line
    for name in record.seq_order:
        yield f"{name} {record.seqs[name]}"
    for name, line in record.gs_lines:
        if name in record.seqs:
            yield line
    for line in record.gc_lines:
        yield line
    for name, line in record.gr_lines:
        if name in record.seqs:
            yield line
    yield "//"


def subset_record(record: StockholmRecord, keep_names: list[str]) -> StockholmRecord:
    keep = set(keep_names)
    return StockholmRecord(
        gf_lines=list(record.gf_lines),
        gc_lines=list(record.gc_lines),
        seq_order=[name for name in record.seq_order if name in keep],
        seqs={name: record.seqs[name] for name in record.seq_order if name in keep},
        gs_lines=[item for item in record.gs_lines if item[0] in keep],
        gr_lines=[item for item in record.gr_lines if item[0] in keep],
        other_lines=list(record.other_lines),
    )


def _append_gc(current: dict, line: str) -> None:
    """Merge ``#=GC`` annotation across interleaved blocks, by tag."""
    parts = line.split(None, 2)
    if len(parts) != 3:
        current["gc_lines"].append(line)
        return
    tag, chunk = parts[1], parts[2]
    for index, existing in enumerate(current["gc_lines"]):
        existing_parts = existing.split(None, 2)
        if len(existing_parts) == 3 and existing_parts[1] == tag:
            current["gc_lines"][index] = f"#=GC {tag} {existing_parts[2]}{chunk}"
            return
    current["gc_lines"].append(line)


def _new_record() -> dict:
    return {
        "gf_lines": [],
        "gc_lines": [],
        "seq_order": [],
        "seqs": {},
        "gs_lines": [],
        "gr_lines": [],
        "other_lines": [],
    }


def _record_from_parts(parts: dict) -> StockholmRecord:
    return StockholmRecord(
        gf_lines=list(parts["gf_lines"]),
        gc_lines=list(parts["gc_lines"]),
        seq_order=list(parts["seq_order"]),
        seqs=dict(parts["seqs"]),
        gs_lines=list(parts["gs_lines"]),
        gr_lines=list(parts["gr_lines"]),
        other_lines=list(parts["other_lines"]),
    )
