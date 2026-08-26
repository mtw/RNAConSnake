#!/usr/bin/env python3

"""Emit ``versions.yaml`` recording every tool version used in a run.

Submission blocker, not a nice-to-have: a calibrated FDR is only reproducible
if the exact toolchain that produced it is recorded alongside it.
"""

from __future__ import annotations

import argparse
import platform
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


# (config tool key, default command, argv used to ask for a version).
# Several of these tools print their version to stderr and/or exit non-zero;
# both are tolerated and the first non-empty output line is kept.
VERSION_PROBES: list[tuple[str, str, list[str]]] = [
    ("rnalalifold", "RNALalifold", ["--version"]),
    ("rnaalifold", "RNAalifold", ["--version"]),
    ("rnafold", "RNAfold", ["--version"]),
    ("rnaz", "RNAz", ["--version"]),
    ("alifoldz", "alifoldz.pl", ["--version"]),
    ("sissiz", "SISSIz", ["--version"]),
    ("rnaz_randomize_aln", "rnazRandomizeAln.pl", ["--version"]),
    ("eslreformat", "esl-reformat", ["-h"]),
    ("rscape", "R-scape", ["-h"]),
    ("cmbuild", "cmbuild", ["-h"]),
    ("cmcalibrate", "cmcalibrate", ["-h"]),
    ("refold", "refold.pl", ["--version"]),
    ("ps2eps", "ps2eps", ["--version"]),
    ("epstopdf", "epstopdf", ["--version"]),
    ("magick", "magick", ["--version"]),
]


VERSION_PATTERN = re.compile(r"\d+\.\d+")


def _first_informative_line(text: str) -> str:
    """Prefer the first line that actually carries a version number.

    Infernal and Easel tools print a banner first and the version on the next
    line, so taking line one alone would report a description instead.
    """
    informative = [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.strip().startswith("Usage")
    ]
    for line in informative[:5]:
        if VERSION_PATTERN.search(line):
            return line
    return informative[0] if informative else ""


def probe(command: str, version_args: list[str], timeout: float = 20.0) -> dict[str, str]:
    tokens = shlex.split(command)
    if not tokens:
        return {"command": command, "path": "MISSING", "version": "unknown"}
    resolved = shutil.which(tokens[0])
    if resolved is None:
        return {"command": command, "path": "MISSING", "version": "unknown"}
    try:
        result = subprocess.run(
            [*tokens, *version_args],
            stdin=subprocess.DEVNULL,
            capture_output=True,
            text=True,
            check=False,
            timeout=timeout,
        )
        output = _first_informative_line((result.stdout or "") + "\n" + (result.stderr or ""))
    except (OSError, subprocess.SubprocessError):
        output = ""
    return {"command": command, "path": resolved, "version": output or "unknown"}


def python_package_version(name: str) -> str:
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover - Python < 3.8 only
        return "unknown"
    try:
        return version(name)
    except PackageNotFoundError:
        return "unknown"


def collect(tools: dict[str, str] | None = None) -> dict[str, object]:
    configured = dict(tools or {})
    external: dict[str, dict[str, str]] = {}
    for key, default, version_args in VERSION_PROBES:
        external[key] = probe(configured.get(key, default), version_args)
    return {
        "rnaconsnake": python_package_version("rnaconsnake"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "snakemake": python_package_version("snakemake"),
        "external_tools": external,
    }


def _yaml_dump(payload: dict[str, object], indent: int = 0) -> str:
    pad = "  " * indent
    lines: list[str] = []
    for key, value in payload.items():
        if isinstance(value, dict):
            lines.append(f"{pad}{key}:")
            lines.append(_yaml_dump(value, indent + 1))
        else:
            text = str(value).replace('"', "'")
            lines.append(f'{pad}{key}: "{text}"')
    return "\n".join(line for line in lines if line)


def write_versions(path: str | Path, tools: dict[str, str] | None = None) -> dict[str, object]:
    payload = collect(tools)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(
        "# RNAcs toolchain versions (generated; do not edit)\n" + _yaml_dump(payload) + "\n",
        encoding="utf-8",
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Record the RNAcs toolchain versions.")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--tool",
        action="append",
        default=[],
        metavar="KEY=COMMAND",
        help="Override a configured tool command. Repeatable.",
    )
    args = parser.parse_args()
    tools = {}
    for item in args.tool:
        key, _, command = item.partition("=")
        if key and command:
            tools[key] = command
    write_versions(args.output, tools)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
