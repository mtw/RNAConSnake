#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

from rnaconsnake.tools.alifold_maxcovar import compute as compute_maxcovar
from rnaconsnake.tools.dereplicate import NR_COLUMNS
from rnaconsnake.workflow_helpers import SUMMARY_FIELDS


def write_json(path: str | Path, payload: dict) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


# The number trailing a label, as RNAz writes it: "Consensus MFE: -18.74".
LABELLED_NUMBER = r"[^0-9-]*(-?\d+(?:\.\d+)?)"

# Labels carrying the consensus MFE, most preferred first. The order is the
# whole point: RNAz prints "Mean single sequence MFE" *before* "Consensus MFE",
# and `re.search` returns the leftmost match, so a single alternation over both
# records the mean single-sequence MFE as the consensus MFE on every real
# output. They are different quantities -- on RNAz 2.1.1 output, -20.93 against
# -18.74 -- so the consensus label is matched on its own, and the bare
# "Mean MFE" is a fallback only for output carrying no consensus line at all.
CONSENSUS_MFE_LABELS = ("consensus MFE", "mean MFE")


def _labelled_number(text: str, labels: tuple[str, ...]) -> str:
    """The number following the first of ``labels`` that ``text`` carries."""
    for label in labels:
        match = re.search(re.escape(label) + LABELLED_NUMBER, text, flags=re.IGNORECASE)
        if match:
            return match.group(1)
    return ""


def cmd_extract_rnaz(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8") if Path(args.input).exists() else ""
    match = re.search(r"probability[^0-9-]*(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    prob = match.group(1) if match else "-1"
    sci_match = re.search(
        r"(?:structure conservation index|SCI)[^0-9-]*(-?\d+(?:\.\d+)?)",
        text,
        flags=re.IGNORECASE,
    )
    write_json(
        args.output,
        {
            "rnazprob": prob,
            "sci": sci_match.group(1) if sci_match else "",
            "consensus_mfe": _labelled_number(text, CONSENSUS_MFE_LABELS),
        },
    )
    return 0


# alifoldz.pl initialises its running minimum to 9999 and prints it unchanged
# when no window was ever scored. It is a sentinel, never a real z-score.
ALIFOLDZ_NO_RESULT = "9999"
ALIFOLDZ_INPUT_LINE = re.compile(r"Input:\s+(\d+)\s+sequences of\s+(\d+)\s+columns")


def cmd_extract_alifoldz(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    score = lines[-1] if lines else ""

    # alifoldz.pl exits 0 even when it parsed nothing at all, so its own report
    # of what it read is the only way to tell a score from a non-result.
    match = ALIFOLDZ_INPUT_LINE.search(text)
    if match and (match.group(1) == "0" or match.group(2) == "0"):
        raise SystemExit(
            f"alifoldz read {match.group(1)} sequences of {match.group(2)} columns from "
            f"{args.input}: it parsed no alignment. alifoldz.pl matches Clustal rows with "
            "[A-Z\\-.]+ and has no /i flag, so a lower-case alignment yields zero sequences."
        )
    if score == ALIFOLDZ_NO_RESULT:
        raise SystemExit(
            f"alifoldz scored no window for {args.input} and emitted its {ALIFOLDZ_NO_RESULT} "
            "sentinel; refusing to record it as a z-score."
        )

    write_json(args.output, {"alifoldzscore": score})
    return 0


def cmd_extract_rscape(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8") if Path(args.input).exists() else ""
    match = re.search(r"^#\s*BPAIRS observed to covary\s+(\d+)\s*$", text, flags=re.MULTILINE)
    count = match.group(1) if match else ""
    write_json(args.output, {"rscape_covary_count": count})
    return 0


# A Clustal sequence line: name, whitespace, then residues/gaps only.
CLUSTAL_SEQUENCE_LINE = re.compile(r"^(\S+)(\s+)([A-Za-z._~-]+)\s*$")


def _uppercase_clustal_sequences(line: str) -> str:
    r"""Upper-case the residues of a Clustal sequence line, leaving the name alone.

    ``refold.pl`` parses alignment rows with ``/^(\S+)\s+([A-Z\-]+)\s*$/`` -- it
    accepts upper case only, as its own "Fixme" admits. Given a lower-case
    alignment (MAFFT output, for instance) it matches nothing, reads zero
    sequences, and then its ``<>`` falls through to STDIN and blocks forever.
    Preparing the alignment for ``refold.pl`` is exactly this step's job, so
    normalise the case here.
    """
    match = CLUSTAL_SEQUENCE_LINE.match(line)
    if not match:
        return line
    name, gap, residues = match.groups()
    return f"{name}{gap}{residues.upper()}"


def cmd_clean_clustal(args: argparse.Namespace) -> int:
    source = Path(args.input).read_text(encoding="utf-8")
    Path(args.backup).write_text(source, encoding="utf-8")

    cleaned_lines = [_uppercase_clustal_sequences(line) for line in source.splitlines() if "*" not in line]
    cleaned = "\n".join(cleaned_lines).rstrip() + "\n"
    Path(args.output).write_text(cleaned, encoding="utf-8")
    return 0


def cmd_extract_consensus(args: argparse.Namespace) -> int:
    """The RNAalifold consensus structure, from its Stockholm output.

    Named for the refold until now, and emitted into the refold directory,
    because it is produced alongside that leg -- but the refold itself is
    `_refold.out`, the per-sequence constrained folds, and nothing here reads
    it. This extracts `#=GC SS_cons` and says so.
    """
    ss_cons = ""
    for line in Path(args.rnaalifold_stk).read_text(encoding="utf-8").splitlines():
        if line.startswith("#=GC SS_cons "):
            parts = line.split(None, 2)
            if len(parts) == 3:
                ss_cons = parts[2]
            break

    write_json(args.output, {"alifold_consstruc": ss_cons})
    return 0


def cmd_run_maxcovar(args: argparse.Namespace) -> int:
    ali_text = Path(args.ali_out).read_text(encoding="utf-8")
    maxcovar, what = compute_maxcovar(ali_text.splitlines())
    Path(args.log).write_text(f"{maxcovar}\t{what}", encoding="utf-8")

    first_ali = ali_text.splitlines()[0].split() if ali_text.splitlines() else []
    nrseq = first_ali[0] if len(first_ali) >= 1 else ""
    alilen = first_ali[5] if len(first_ali) >= 6 else ""

    write_json(
        args.output,
        {
            "maxcovarval": str(maxcovar),
            "nrseq": nrseq,
            "alilen": alilen,
            "maxcovarcount": what,
        },
    )
    return 0


def cmd_write_disabled(args: argparse.Namespace) -> int:
    write_json(args.output, {args.key: args.value})
    return 0


def cmd_combine_summary(args: argparse.Namespace) -> int:
    payload: dict[str, str] = {"wbn": args.wbn}
    for src in args.inputs:
        payload.update(json.loads(Path(src).read_text(encoding="utf-8")))
    write_json(args.output, payload)
    return 0


def _load_summary_records(paths: list[str]) -> list[dict[str, str]]:
    records = [json.loads(Path(path).read_text(encoding="utf-8")) for path in sorted(paths)]
    normalized: list[dict[str, str]] = []
    for record in records:
        normalized.append({field: str(record.get(field, "")) for field in SUMMARY_FIELDS})
    return normalized


def _numeric_value(value: str, default: float = float("-inf")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _report_sort_key(record: dict[str, str]) -> tuple[float, float]:
    """Report order: strongest covariation first, then strongest AlifoldZ.

    AlifoldZ is negated because a *more negative* z-score is the more
    significant one -- the same ordering de-replication uses to pick a locus
    representative. Missing values rank last in both positions.
    """
    return (
        _numeric_value(record.get("maxcovarval", "")),
        -_numeric_value(record.get("alifoldzscore", ""), default=float("inf")),
    )


def cmd_render_reports(args: argparse.Namespace) -> int:
    records = _load_summary_records(args.inputs)
    records = sorted(records, key=_report_sort_key, reverse=True)
    for path in [args.log, args.csv, args.markdown]:
        if path:
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    with open(args.log, "w", encoding="utf-8") as log_handle:
        log_handle.writelines(
            "X {wbn} {nrseq} {alilen} maxcovar {maxcovarval} {maxcovarcount} rscape "
            "{rscape_covary_count} {rnazprob} {sci} {consensus_mfe} {alifoldzscore} "
            "{alifold_consstruc}\n".format(**record)
            for record in records
        )

    with open(args.csv, "w", encoding="utf-8", newline="") as csv_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in SUMMARY_FIELDS})

    if args.markdown:
        _write_markdown_table(
            args.markdown,
            heading=f"# RNAConSnake Summary: {args.label}",
            preamble=[f"- Records: {len(records)}"],
            columns=SUMMARY_FIELDS,
            rows=records,
        )
    return 0


def _markdown_cell(value: object) -> str:
    return str(value if value is not None else "").replace("\n", " ").replace("|", "\\|")


def _markdown_table(columns: list[str], rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(_markdown_cell(row.get(c, "")) for c in columns) + " |")
    return lines


def _write_markdown_table(path, heading, preamble, columns, rows) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(heading + "\n\n")
        handle.writelines(line + "\n" for line in preamble)
        handle.write("\n")
        handle.writelines(line + "\n" for line in _markdown_table(columns, rows))


def _read_csv_rows(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def cmd_render_markdown(args: argparse.Namespace) -> int:
    """Render the human-readable summary: non-redundant loci first, then every window.

    The de-replicated block comes first because that is the list to read: one
    row per locus rather than one per RNALalifold window, which reports a single
    element repeatedly (see docs/dereplication.md). The full per-window table
    follows unchanged, so nothing is hidden.
    """
    nr_rows = _read_csv_rows(args.nr)
    full_rows = _read_csv_rows(args.full)

    # Members last: it is the widest column and the least often read.
    nr_columns = [c for c in NR_COLUMNS if c != "members"] + ["members"]

    lines: list[str] = [f"# RNAConSnake Summary: {args.label}", ""]
    lines.append(f"- Loci (non-redundant): {len(nr_rows)}")
    lines.append(f"- Windows (all): {len(full_rows)}")
    if args.method:
        lines.append(f"- De-replication method: `{args.method}`")
    lines += [
        "",
        "## Non-redundant candidates",
        "",
        "One row per locus: the best-scoring window, and the windows collapsed "
        "into it. `locus_start`/`locus_end` is the union of the members, not the "
        "extent of the reported element. See `docs/dereplication.md`.",
        "",
    ]
    lines += _markdown_table(nr_columns, nr_rows)
    lines += [
        "",
        "## All windows",
        "",
        f"Every RNALalifold window, including the {max(len(full_rows) - len(nr_rows), 0)} collapsed above.",
        "",
    ]
    lines += _markdown_table(SUMMARY_FIELDS, full_rows)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RNAConSnake post-processing helpers.")
    sub = parser.add_subparsers(dest="command", required=True)

    rnaz = sub.add_parser("extract-rnaz")
    rnaz.add_argument("--input", required=True)
    rnaz.add_argument("--output", required=True)
    rnaz.set_defaults(func=cmd_extract_rnaz)

    alifoldz = sub.add_parser("extract-alifoldz")
    alifoldz.add_argument("--input", required=True)
    alifoldz.add_argument("--output", required=True)
    alifoldz.set_defaults(func=cmd_extract_alifoldz)

    rscape = sub.add_parser("extract-rscape")
    rscape.add_argument("--input", required=True)
    rscape.add_argument("--output", required=True)
    rscape.set_defaults(func=cmd_extract_rscape)

    clean = sub.add_parser("clean-clustal")
    clean.add_argument("--input", required=True)
    clean.add_argument("--backup", required=True)
    clean.add_argument("--output", required=True)
    clean.set_defaults(func=cmd_clean_clustal)

    consensus = sub.add_parser("extract-consensus")
    consensus.add_argument("--rnaalifold-stk", required=True)
    consensus.add_argument("--output", required=True)
    consensus.set_defaults(func=cmd_extract_consensus)

    maxcov = sub.add_parser("run-maxcovar")
    maxcov.add_argument("--ali-out", required=True)
    maxcov.add_argument("--log", required=True)
    maxcov.add_argument("--output", required=True)
    maxcov.set_defaults(func=cmd_run_maxcovar)

    disabled = sub.add_parser("write-disabled")
    disabled.add_argument("--key", required=True)
    disabled.add_argument("--value", required=True)
    disabled.add_argument("--output", required=True)
    disabled.set_defaults(func=cmd_write_disabled)

    combine = sub.add_parser("combine-summary")
    combine.add_argument("--wbn", required=True)
    combine.add_argument("--output", required=True)
    combine.add_argument("inputs", nargs="+")
    combine.set_defaults(func=cmd_combine_summary)

    reports = sub.add_parser("write-summary-outputs")
    reports.add_argument("--label", required=True)
    reports.add_argument("--log", required=True)
    reports.add_argument("--csv", required=True)
    reports.add_argument("--markdown")
    reports.add_argument("inputs", nargs="+")
    reports.set_defaults(func=cmd_render_reports)

    markdown = sub.add_parser("render-markdown")
    markdown.add_argument("--label", required=True)
    markdown.add_argument("--nr", required=True, help="Non-redundant locus CSV.")
    markdown.add_argument("--full", required=True, help="Per-window summary CSV.")
    markdown.add_argument("--output", required=True)
    markdown.add_argument("--method", help="De-replication method, for the header.")
    markdown.set_defaults(func=cmd_render_markdown)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
