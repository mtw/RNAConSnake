#!/usr/bin/env python3

"""Fold and plot an arbitrary alignment region.

The screen reports the *best-scoring* window of a locus, which can be much
shorter than the element it belongs to -- a dumbbell whose 5' lobe scores best
is reported as that lobe. For curation you need the whole element folded and
drawn, not the fragment the screen preferred.

This extracts a column range from an alignment, puts it through the same
preprocessing and the same RNAalifold invocation the workflow uses, and emits
the consensus structure, the plots and the scores for that exact span.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from rnaconsnake.tools.alignment_io import Alignment, read_stockholm_alignment
from rnaconsnake.workflow_helpers import normalize_rnaalifold_side_output, run_checked


def slice_alignment(alignment: Alignment, start: int, end: int) -> Alignment:
    """1-based inclusive column slice."""
    if start < 1 or end > alignment.length or start > end:
        raise ValueError(f"span {start}-{end} is outside the alignment (1-{alignment.length})")
    lo, hi = start - 1, end
    return Alignment(
        order=list(alignment.order),
        seqs={name: alignment.seqs[name][lo:hi] for name in alignment.order},
    )


def write_stockholm(alignment: Alignment, path: Path, identifier: str) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        handle.write(f"#=GF ID {identifier}\n")
        handle.writelines(f"{name} {alignment.seqs[name]}\n" for name in alignment.order)
        handle.write("//\n")


def preprocess(source: Path, workdir: Path, label: str, gapratio: float, max_n: int) -> Path:
    """The workflow's own cleanup: drop gap-dominated and redundant sequences."""
    remgap = workdir / f"{label}.remgap.stk"
    run_checked(
        [
            "python3",
            "-m",
            "rnaconsnake.tools.remove_gaponly",
            "-a",
            str(source),
            "-i",
            "stockholm",
            "-r",
            str(gapratio),
            "-n",
            str(max_n),
        ],
        stdout_path=remgap,
    )
    cleaned = workdir / f"{label}.stk"
    run_checked(
        ["python3", "-m", "rnaconsnake.tools.strip_aln", "-a", str(remgap), "-f", "S", "--nosingle"],
        stdout_path=cleaned,
    )
    return cleaned


def fold_and_plot(cleaned: Path, outdir: Path, label: str, tools: dict[str, str]) -> Path:
    """RNAalifold with the workflow's parameters, then EPS and PDF."""
    run_checked(
        [
            *shlex.split(tools["rnaalifold"]),
            "-t4",
            "--aln",
            "--color",
            "-r",
            "--cfactor",
            "0.6",
            "--nfactor",
            "0.5",
            "-p",
            "--aln-EPS-cols=200",
            "--aln-stk",
            "-f",
            "S",
            "--id-prefix",
            label,
        ],
        stdin_path=cleaned,
        stdout_path=outdir / f"{label}.alifold.out",
        cwd=str(outdir),
    )
    for stray in [outdir / "alirna.ps"]:
        if stray.exists():
            stray.unlink()
    default_stk = outdir / "RNAalifold_results.stk"
    result_stk = outdir / f"{label}.RNAalifold_results.stk"
    if default_stk.exists() and not result_stk.exists():
        default_stk.rename(result_stk)
    for suffix in ["_ali.out", "_dp.ps", "_aln.ps", "_ss.ps"]:
        normalize_rnaalifold_side_output(outdir, outdir / f"{label}{suffix}", suffix)
    for base in ["_aln", "_ss"]:
        run_checked([*shlex.split(tools["ps2eps"]), f"{label}{base}.ps"], cwd=str(outdir))
        run_checked([*shlex.split(tools["epstopdf"]), f"{label}{base}.eps"], cwd=str(outdir))
    return result_stk


def score(cleaned: Path, outdir: Path, label: str, tools: dict[str, str]) -> dict:
    """RNAz and AlifoldZ on this exact span, for comparison with the screen."""
    clustal = outdir / f"{label}.aln"
    run_checked([*shlex.split(tools["eslreformat"]), "clustal", str(cleaned)], stdout_path=clustal)

    rnaz_txt = outdir / f"{label}.rnaz.txt"
    run_checked([*shlex.split(tools["rnaz"]), "-d", "-n", str(clustal)], stdout_path=rnaz_txt)
    rnaz_json = outdir / f"{label}.rnaz.json"
    run_checked(
        [
            "python3",
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "extract-rnaz",
            "--input",
            str(rnaz_txt),
            "--output",
            str(rnaz_json),
        ]
    )
    scores = json.loads(rnaz_json.read_text(encoding="utf-8"))

    with open(clustal, encoding="utf-8") as handle:
        # In outdir, not wherever the user invoked this: alifoldz.pl shells out
        # to RNAalifold without --noPS and drops an undeclared "alirna.ps" into
        # the directory it runs in.
        result = subprocess.run(
            [*shlex.split(tools["alifoldz"]), "-f", "-t", "0.0"],
            stdin=handle,
            capture_output=True,
            text=True,
            check=False,
            cwd=str(outdir),
        )
    stray = outdir / "alirna.ps"
    if stray.exists():
        stray.unlink()
    lines = [line.strip() for line in (result.stdout or "").splitlines() if line.strip()]
    value = lines[-1] if lines else ""
    scores["alifoldzscore"] = "NA" if (result.returncode != 0 or value == "9999") else value
    return scores


def main() -> int:
    parser = argparse.ArgumentParser(description="Fold and plot one alignment region.")
    parser.add_argument("--alignment", required=True)
    parser.add_argument("--start", type=int, required=True, help="1-based inclusive column")
    parser.add_argument("--end", type=int, required=True, help="1-based inclusive column")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--remove-gaponly-gapratio", type=float, default=0.5)
    parser.add_argument("--remove-gaponly-max-n", type=int, default=0)
    parser.add_argument("--rnaalifold", default="RNAalifold")
    parser.add_argument("--rnaz", default="RNAz")
    parser.add_argument("--alifoldz", default="alifoldz.pl")
    parser.add_argument("--eslreformat", default="esl-reformat")
    parser.add_argument("--ps2eps", default="ps2eps")
    parser.add_argument("--epstopdf", default="epstopdf")
    parser.add_argument("--no-scores", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    tools = {
        "rnaalifold": args.rnaalifold,
        "rnaz": args.rnaz,
        "alifoldz": args.alifoldz,
        "eslreformat": args.eslreformat,
        "ps2eps": args.ps2eps,
        "epstopdf": args.epstopdf,
    }

    alignment = read_stockholm_alignment(args.alignment)
    region = slice_alignment(alignment, args.start, args.end)
    raw = outdir / f"{args.label}.region.stk"
    write_stockholm(region, raw, args.label)
    cleaned = preprocess(raw, outdir, args.label, args.remove_gaponly_gapratio, args.remove_gaponly_max_n)
    fold_and_plot(cleaned, outdir, args.label, tools)

    payload = {
        "label": args.label,
        "start": args.start,
        "end": args.end,
        "columns": args.end - args.start + 1,
        "n_seq_after_cleanup": read_stockholm_alignment(cleaned).n_seq,
    }
    if not args.no_scores:
        payload.update(score(cleaned, outdir, args.label, tools))
    (outdir / f"{args.label}.scores.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"{args.label}: {args.start}-{args.end} ({payload['columns']} cols) -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
