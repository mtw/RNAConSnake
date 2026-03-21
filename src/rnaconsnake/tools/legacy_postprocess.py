#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import sys
from pathlib import Path

from rnaconsnake.tools.alifold_maxcovar import compute as compute_maxcovar


def write_json(path: str | Path, payload: dict) -> None:
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def cmd_extract_rnaz(args: argparse.Namespace) -> int:
    text = Path(args.input).read_text(encoding="utf-8") if Path(args.input).exists() else ""
    match = re.search(r"probability[^0-9-]*(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    prob = match.group(1) if match else "-1"
    write_json(args.output, {"rnazprob": prob})
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
    "alifoldzscore",
    "refold_firstseq",
    "refold_firststruc",
    "alifold_consstruc",
]

HTML_FIELDS = [
    field
    for field in SUMMARY_FIELDS
    if field not in {"refold_firstseq", "refold_firststruc"}
]

NUMERIC_SORT_FIELDS = {
    "nrseq",
    "alilen",
    "maxcovarval",
    "maxcovarcount",
    "rscape_covary_count",
    "rnazprob",
    "alifoldzscore",
}


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


def _html_cell_content(field: str, record: dict[str, str], html_path: str, label: str) -> str:
    value = str(record.get(field, "")).replace("\n", " ")
    escaped_value = html.escape(value)
    if field == "rscape_covary_count":
        wbn = str(record.get("wbn", "")).replace("\n", " ")
        html_dir = Path(html_path).parent
        sto_pdf = html_dir.parent.parent / "rscape" / label / f"{wbn}.sto.pdf"
        sto_href = html.escape(os.path.relpath(sto_pdf, html_dir)) if sto_pdf.exists() else ""
        if sto_href and value not in {"", "0"}:
            return f"<a class='rscape-link' href='{sto_href}' data-preview-src='{sto_href}'><code>{escaped_value}</code></a>"
        return f"<code>{escaped_value}</code>"

    if field != "wbn":
        return f"<code>{escaped_value}</code>"

    html_dir = Path(html_path).parent
    rnaalifold_dir = html_dir.parent.parent / "rnaalifold" / label / value
    aln_pdf = rnaalifold_dir / f"{value}_aln.pdf"
    ss_pdf = rnaalifold_dir / f"{value}_ss.pdf"
    aln_href = html.escape(os.path.relpath(aln_pdf, html_dir)) if aln_pdf.exists() else ""
    ss_href = html.escape(os.path.relpath(ss_pdf, html_dir)) if ss_pdf.exists() else ""
    preview_attr = f" data-preview-src='{aln_href}'" if aln_href else ""
    links = []
    if aln_href:
        links.append(f"<a class='plot-link' href='{aln_href}'>aln pdf</a>")
    if ss_href:
        links.append(f"<a class='plot-link' href='{ss_href}'>ss pdf</a>")
    link_html = f"<span class='plot-links'>{' '.join(links)}</span>" if links else ""
    if aln_href:
        label_html = f"<a class='wbn-link' href='{aln_href}'{preview_attr}>{escaped_value}</a>"
    else:
        label_html = f"<span class='wbn-link'>{escaped_value}</span>"
    return f"<span class='wbn-cell'>{label_html}{link_html}</span>"


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
    for path in [args.log, args.csv, args.markdown, args.html]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)

    with open(args.log, "w", encoding="utf-8") as log_handle:
        for record in records:
            log_handle.write(
                "X {wbn} {nrseq} {alilen} maxcovar {maxcovarval} {maxcovarcount} rscape {rscape_covary_count} {rnazprob} {alifoldzscore} {refold_firstseq} {refold_firststruc} {alifold_consstruc}\n".format(
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

    with open(args.html, "w", encoding="utf-8") as html_handle:
        html_handle.write(
            "<!doctype html>\n"
            "<html lang='en'>\n"
            "<head>\n"
            "  <meta charset='utf-8'>\n"
            f"  <title>RNAConSnake Summary: {html.escape(args.label)}</title>\n"
            "  <style>\n"
            "    body { font-family: Helvetica, Arial, sans-serif; margin: 2rem; }\n"
            "    h1 { margin-bottom: 0.25rem; }\n"
            "    .meta { color: #555; margin-bottom: 1rem; }\n"
            "    .controls { display: grid; grid-template-columns: repeat(4, minmax(12rem, 1fr)); gap: 0.75rem; margin-bottom: 1rem; }\n"
            "    .control { display: flex; flex-direction: column; gap: 0.25rem; }\n"
            "    .control label { font-size: 0.85rem; color: #444; }\n"
            "    .control input { padding: 0.35rem 0.5rem; border: 1px solid #c8c8c8; border-radius: 4px; }\n"
            "    table { border-collapse: collapse; width: 100%; font-size: 0.9rem; }\n"
            "    th, td { border: 1px solid #d0d0d0; padding: 0.35rem 0.5rem; text-align: left; vertical-align: top; }\n"
            "    thead th { position: sticky; top: 0; background: #f5f5f5; }\n"
            "    th button { all: unset; cursor: pointer; font-weight: 700; }\n"
            "    tbody tr:nth-child(even) { background: #fafafa; }\n"
            "    code { font-family: Menlo, Consolas, monospace; font-size: 0.85em; }\n"
            "    .wbn-cell { display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap; }\n"
            "    .wbn-link { color: #12436d; text-decoration: underline; cursor: pointer; }\n"
            "    .plot-links { display: inline-flex; gap: 0.35rem; }\n"
            "    .plot-link { font-size: 0.8rem; color: #555; text-decoration: none; border: 1px solid #c8c8c8; border-radius: 999px; padding: 0.05rem 0.4rem; }\n"
            "    .plot-link:hover { background: #efefef; }\n"
            "    .preview { position: fixed; right: 1rem; bottom: 1rem; width: min(42vw, 560px); height: min(60vh, 720px); background: rgba(255,255,255,0.98); border: 1px solid #d0d0d0; box-shadow: 0 8px 24px rgba(0,0,0,0.12); padding: 0.5rem; display: none; }\n"
            "    .preview iframe { display: block; width: 100%; height: 100%; border: 0; }\n"
            "  </style>\n"
            "</head>\n"
            "<body>\n"
            f"  <h1>RNAConSnake Summary: {html.escape(args.label)}</h1>\n"
            "  <div class='meta'>"
            f"Records: {len(records)}. "
            "Default sort: maxcovarval desc, then alifoldzscore desc."
            "</div>\n"
            "  <div class='controls'>\n"
            "    <div class='control'><label for='filter-wbn'>wbn contains</label><input id='filter-wbn' type='search' placeholder='RC_100_0001...'></div>\n"
            "    <div class='control'><label for='filter-nrseq'>nrseq minimum (>=)</label><input id='filter-nrseq' type='number' step='any'></div>\n"
            "    <div class='control'><label for='filter-alilen'>alilen minimum (>=)</label><input id='filter-alilen' type='number' step='any' value='20'></div>\n"
            "    <div class='control'><label for='filter-maxcovarval'>maxcovarval minimum (>=)</label><input id='filter-maxcovarval' type='number' step='any'></div>\n"
            "    <div class='control'><label for='filter-maxcovarcount'>maxcovarcount minimum (>=)</label><input id='filter-maxcovarcount' type='number' step='any'></div>\n"
            "    <div class='control'><label for='filter-rscape_covary_count'>rscape covary minimum (>=)</label><input id='filter-rscape_covary_count' type='number' step='any'></div>\n"
            "    <div class='control'><label for='filter-rnazprob'>rnazprob minimum (>=)</label><input id='filter-rnazprob' type='number' step='any' value='0.9'></div>\n"
            "    <div class='control'><label for='filter-alifoldzscore'>alifoldzscore maximum (&lt;=)</label><input id='filter-alifoldzscore' type='number' step='any' value='-1'></div>\n"
            "  </div>\n"
            "  <div id='preview' class='preview'><iframe id='preview-frame' title='alignment preview'></iframe></div>\n"
            "  <table>\n"
            "    <thead><tr>"
            + "".join(
                (
                    f"<th data-field='{html.escape(field)}' data-sort-type='"
                    f"{'number' if field in NUMERIC_SORT_FIELDS else 'text'}'>"
                    f"<button type='button'>{html.escape(field)}</button></th>"
                )
                for field in HTML_FIELDS
            )
            + "</tr></thead>\n"
            "    <tbody>\n"
        )
        for record in records:
            html_handle.write(
                "      <tr>"
                + "".join(
                    (
                        f"<td data-value='{html.escape(str(record.get(field, '')).replace(chr(10), ' '))}'>"
                        f"{_html_cell_content(field, record, args.html, args.label)}</td>"
                    )
                    for field in HTML_FIELDS
                )
                + "</tr>\n"
            )
        html_handle.write(
            "    </tbody>\n"
            "  </table>\n"
            "  <script>\n"
            "    const table = document.querySelector('table');\n"
            "    const tbody = table.querySelector('tbody');\n"
            "    const headers = Array.from(table.querySelectorAll('thead th'));\n"
            "    const preview = document.getElementById('preview');\n"
            "    const previewFrame = document.getElementById('preview-frame');\n"
            "    const filterIds = {\n"
            "      wbn: document.getElementById('filter-wbn'),\n"
            "      nrseq: document.getElementById('filter-nrseq'),\n"
            "      alilen: document.getElementById('filter-alilen'),\n"
            "      maxcovarval: document.getElementById('filter-maxcovarval'),\n"
            "      maxcovarcount: document.getElementById('filter-maxcovarcount'),\n"
            "      rscape_covary_count: document.getElementById('filter-rscape_covary_count'),\n"
            "      rnazprob: document.getElementById('filter-rnazprob'),\n"
            "      alifoldzscore: document.getElementById('filter-alifoldzscore')\n"
            "    };\n"
            "    const rows = () => Array.from(tbody.querySelectorAll('tr'));\n"
            "    const fieldIndex = Object.fromEntries(headers.map((header, index) => [header.dataset.field, index]));\n"
            "    const compareValues = (a, b, type, direction) => {\n"
            "      if (type === 'number') {\n"
            "        const av = Number.parseFloat(a);\n"
            "        const bv = Number.parseFloat(b);\n"
            "        const an = Number.isNaN(av) ? Number.NEGATIVE_INFINITY : av;\n"
            "        const bn = Number.isNaN(bv) ? Number.NEGATIVE_INFINITY : bv;\n"
            "        return direction * (an - bn);\n"
            "      }\n"
            "      return direction * a.localeCompare(b);\n"
            "    };\n"
            "    const sortBy = (specs) => {\n"
            "      const sorted = rows().sort((rowA, rowB) => {\n"
            "        for (const spec of specs) {\n"
            "          const cellA = rowA.children[spec.index].dataset.value || '';\n"
            "          const cellB = rowB.children[spec.index].dataset.value || '';\n"
            "          const delta = compareValues(cellA, cellB, spec.type, spec.direction);\n"
            "          if (delta !== 0) return delta;\n"
            "        }\n"
            "        return 0;\n"
            "      });\n"
            "      tbody.replaceChildren(...sorted);\n"
            "    };\n"
            "    const passesFilters = (row) => {\n"
            "      const wbnNeedle = (filterIds.wbn.value || '').toLowerCase().trim();\n"
            "      if (wbnNeedle) {\n"
            "        const wbn = (row.children[fieldIndex.wbn].dataset.value || '').toLowerCase();\n"
            "        if (!wbn.includes(wbnNeedle)) return false;\n"
            "      }\n"
            "      for (const field of ['nrseq', 'alilen', 'maxcovarval', 'maxcovarcount', 'rscape_covary_count', 'rnazprob', 'alifoldzscore']) {\n"
            "        const needle = filterIds[field].value;\n"
            "        if (!needle) continue;\n"
            "        const cell = row.children[fieldIndex[field]].dataset.value || '';\n"
            "        const actual = Number.parseFloat(cell);\n"
            "        const threshold = Number.parseFloat(needle);\n"
            "        if (Number.isNaN(actual)) return false;\n"
            "        if (field === 'alifoldzscore') {\n"
            "          if (actual > threshold) return false;\n"
            "        } else if (actual < threshold) {\n"
            "          return false;\n"
            "        }\n"
            "      }\n"
            "      return true;\n"
            "    };\n"
            "    const applyFilters = () => {\n"
            "      rows().forEach((row) => {\n"
            "        row.style.display = passesFilters(row) ? '' : 'none';\n"
            "      });\n"
            "    };\n"
            "    headers.forEach((header, index) => {\n"
            "      let descending = true;\n"
            "      header.querySelector('button').addEventListener('click', () => {\n"
            "        sortBy([{ index, type: header.dataset.sortType, direction: descending ? -1 : 1 }]);\n"
            "        applyFilters();\n"
            "        descending = !descending;\n"
            "      });\n"
            "    });\n"
            "    table.querySelectorAll('[data-preview-src]').forEach((link) => {\n"
            "      link.addEventListener('mouseenter', () => {\n"
            "        previewFrame.src = link.dataset.previewSrc;\n"
            "        preview.style.display = 'block';\n"
            "      });\n"
            "      link.addEventListener('mouseleave', () => {\n"
            "        preview.style.display = 'none';\n"
            "        previewFrame.removeAttribute('src');\n"
            "      });\n"
            "    });\n"
            "    Object.values(filterIds).forEach((input) => input.addEventListener('input', applyFilters));\n"
            "    sortBy([\n"
            "      { index: fieldIndex.maxcovarval, type: 'number', direction: -1 },\n"
            "      { index: fieldIndex.alifoldzscore, type: 'number', direction: -1 }\n"
            "    ]);\n"
            "    applyFilters();\n"
            "  </script>\n"
            "</body>\n"
            "</html>\n"
        )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Legacy pp_RNALalifold post-processing helpers.")
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

    reports = sub.add_parser("render-reports")
    reports.add_argument("--label", required=True)
    reports.add_argument("--log", required=True)
    reports.add_argument("--csv", required=True)
    reports.add_argument("--markdown", required=True)
    reports.add_argument("--html", required=True)
    reports.add_argument("inputs", nargs="+")
    reports.set_defaults(func=cmd_render_reports)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
