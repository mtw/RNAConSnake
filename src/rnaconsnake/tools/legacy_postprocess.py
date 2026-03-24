#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

from rnaconsnake.tools.alifold_maxcovar import compute as compute_maxcovar


def write_json(path: str | Path, payload: dict) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def cmd_extract_rnaz(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8") if Path(args.input).exists() else ""
    match = re.search(r"probability[^0-9-]*(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    prob = match.group(1) if match else "-1"
    sci_match = re.search(
        r"(?:structure conservation index|SCI)[^0-9-]*(-?\d+(?:\.\d+)?)",
        text,
        flags=re.IGNORECASE,
    )
    mfe_match = re.search(
        r"(?:consensus MFE|mean single sequence MFE|mean MFE)[^0-9-]*(-?\d+(?:\.\d+)?)",
        text,
        flags=re.IGNORECASE,
    )
    write_json(
        args.output,
        {
            "rnazprob": prob,
            "sci": sci_match.group(1) if sci_match else "",
            "consensus_mfe": mfe_match.group(1) if mfe_match else "",
        },
    )
    return 0


def cmd_extract_alifoldz(args: argparse.Namespace) -> int:
    lines = [line.strip() for line in Path(args.input).read_text(encoding="utf-8").splitlines() if line.strip()]
    score = lines[-1] if lines else ""
    write_json(args.output, {"alifoldzscore": score})
    return 0


def cmd_extract_rscape(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8") if Path(args.input).exists() else ""
    match = re.search(r"^#\s*BPAIRS observed to covary\s+(\d+)\s*$", text, flags=re.MULTILINE)
    count = match.group(1) if match else ""
    write_json(args.output, {"rscape_covary_count": count})
    return 0


def cmd_clean_clustal(args: argparse.Namespace) -> int:
    source = Path(args.input).read_text(encoding="utf-8")
    Path(args.backup).write_text(source, encoding="utf-8")

    cleaned_lines = [line for line in source.splitlines() if "*" not in line]
    cleaned = "\n".join(cleaned_lines).rstrip() + "\n"
    Path(args.output).write_text(cleaned, encoding="utf-8")
    return 0


def cmd_extract_refold(args: argparse.Namespace) -> int:
    refold_lines = Path(args.refold_output).read_text(encoding="utf-8").splitlines()
    firstseq = refold_lines[1] if len(refold_lines) >= 2 else ""
    firststruc = refold_lines[2] if len(refold_lines) >= 3 else ""

    ss_cons = ""
    for line in Path(args.rnaalifold_stk).read_text(encoding="utf-8").splitlines():
        if line.startswith("#=GC SS_cons "):
            parts = line.split(None, 2)
            if len(parts) == 3:
                ss_cons = parts[2]
            break

    write_json(
        args.output,
        {
            "refold_firstseq": firstseq,
            "refold_firststruc": firststruc,
            "alifold_consstruc": ss_cons,
        },
    )
    return 0


def cmd_run_maxcovar(args: argparse.Namespace) -> int:
    ali_text = Path(args.ali_out).read_text(encoding="utf-8")
    maxcovar, what = compute_maxcovar(ali_text.splitlines())
    Path(args.log).write_text(f"{maxcovar}\t{what}", encoding="utf-8")

    first_ali = ali_text.splitlines()[0].split() if ali_text.splitlines() else []
    nrseq = first_ali[0] if len(first_ali) >= 1 else ""
    alilen = first_ali[5] if len(first_ali) >= 6 else ""

    first_log = [str(maxcovar), what]
    maxcovarcount = first_log[1] if len(first_log) >= 2 else ""

    write_json(
        args.output,
        {
            "maxcovarval": str(maxcovar),
            "nrseq": nrseq,
            "alilen": alilen,
            "maxcovarcount": maxcovarcount,
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


SUMMARY_FIELDS = [
    "wbn",
    "nrseq",
    "alilen",
    "maxcovarval",
    "maxcovarcount",
    "rscape_covary_count",
    "rnazprob",
    "sci",
    "consensus_mfe",
    "alifoldzscore",
    "alifold_consstruc",
]

def _load_summary_records(paths: list[str]) -> list[dict[str, str]]:
    records = [json.loads(Path(path).read_text(encoding="utf-8")) for path in sorted(paths)]
    normalized: list[dict[str, str]] = []
    for record in records:
        normalized.append({field: str(record.get(field, "")) for field in SUMMARY_FIELDS})
    return normalized


def _numeric_value(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("-inf")


def cmd_render_reports(args: argparse.Namespace) -> int:
    records = _load_summary_records(args.inputs)
    records = sorted(
        records,
        key=lambda record: (
            _numeric_value(record.get("maxcovarval", "")),
            _numeric_value(record.get("alifoldzscore", "")),
        ),
        reverse=True,
    )
    for path in [args.log, args.csv, args.markdown]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    with open(args.log, "w", encoding="utf-8") as log_handle:
        for record in records:
            log_handle.write(
                "X {wbn} {nrseq} {alilen} maxcovar {maxcovarval} {maxcovarcount} rscape {rscape_covary_count} {rnazprob} {sci} {consensus_mfe} {alifoldzscore} {alifold_consstruc}\n".format(
                    **record
                )
            )

    with open(args.csv, "w", encoding="utf-8", newline="") as csv_handle:
        writer = csv.DictWriter(csv_handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record.get(field, "") for field in SUMMARY_FIELDS})

    with open(args.markdown, "w", encoding="utf-8") as md_handle:
        md_handle.write(f"# RNAConSnake Summary: {args.label}\n\n")
        md_handle.write(f"- Records: {len(records)}\n\n")
        md_handle.write("| " + " | ".join(SUMMARY_FIELDS) + " |\n")
        md_handle.write("| " + " | ".join(["---"] * len(SUMMARY_FIELDS)) + " |\n")
        for record in records:
            row = [str(record.get(field, "")).replace("\n", " ") for field in SUMMARY_FIELDS]
            md_handle.write("| " + " | ".join(row) + " |\n")
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

    refold = sub.add_parser("extract-refold")
    refold.add_argument("--refold-output", required=True)
    refold.add_argument("--rnaalifold-stk", required=True)
    refold.add_argument("--output", required=True)
    refold.set_defaults(func=cmd_extract_refold)

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
    reports.add_argument("--markdown", required=True)
    reports.add_argument("inputs", nargs="+")
    reports.set_defaults(func=cmd_render_reports)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
