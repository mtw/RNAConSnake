from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from importlib.resources import files
from pathlib import Path

import yaml

from rnaconsnake import __version__
from rnaconsnake.export_bundle import build_export

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback only
    tqdm = None

DEFAULT_RUNTIME_DEPENDENCIES = [
    "RNALalifold",
    "snakemake",
    "esl-reformat",
    "RNAz",
    "alifoldz.pl",
    "RNAalifold",
    "ps2eps",
    "epstopdf",
    "refold.pl",
    "RNAfold",
    "magick",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the packaged RNAConSnake Snakemake workflow.")
    parser.add_argument("--cores", default="all", help="Snakemake cores value, e.g. 4 or all")
    parser.add_argument("--latency-wait", default="20", help="Snakemake latency wait")
    parser.add_argument("--snakefile", help="Override packaged snakefile path")
    parser.add_argument(
        "--maxbpspan",
        type=int,
        action="append",
        help="Override maxbpspan values from config.yaml. Repeat the option to provide multiple window sizes.",
    )
    parser.add_argument(
        "--rscape",
        action="store_true",
        help="Enable optional per-hit R-scape analysis.",
    )
    parser.add_argument(
        "--rnaz-shuffle",
        action="store_true",
        help="Enable RNAz alignment shuffling by omitting the default -n flag.",
    )
    parser.add_argument(
        "--output-dir",
        help="Working directory where all workflow outputs will be written.",
    )
    parser.add_argument(
        "--input-alignment",
        help="Run the workflow on an arbitrary input Stockholm alignment file instead of data/<sample>.stk.",
    )
    parser.add_argument(
        "--export-bundle",
        help="Optional output directory for a post-run RNAConSnake export bundle.",
    )
    parser.add_argument(
        "--export-dataset-id",
        help="Optional dataset_id override for the export bundle.",
    )
    parser.add_argument(
        "--export-dataset-label",
        help="Optional dataset_label override for the export bundle.",
    )
    parser.add_argument(
        "--export-feature-id",
        help="Optional feature_id override for the export bundle.",
    )
    parser.add_argument(
        "--export-feature-label",
        help="Optional feature_label override for the export bundle.",
    )
    parser.add_argument(
        "--export-feature-type",
        default="other",
        help="Feature type for the export bundle, default: other.",
    )
    parser.add_argument(
        "--export-source-label",
        help="Optional source/provenance label for the export bundle.",
    )
    parser.add_argument(
        "--export-description",
        help="Optional dataset/feature description for the export bundle.",
    )
    parser.add_argument(
        "--conservative",
        action="store_true",
        help="Trust existing outputs more conservatively by only rerunning on mtime-based staleness.",
    )
    parser.add_argument(
        "--progress",
        action="store_true",
        help="Render interactive tqdm progress bars while Snakemake runs.",
    )
    parser.add_argument(
        "--no-progress",
        action="store_true",
        help="Disable interactive tqdm progress bars and forward Snakemake output directly.",
    )
    parser.add_argument(
        "--show-tool-paths",
        action="store_true",
        help="Show the configured external tool commands and resolved executable paths, then exit.",
    )
    parser.add_argument("--check-deps", action="store_true", help="Check external runtime dependencies and exit")
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    parser.add_argument("snakemake_args", nargs=argparse.REMAINDER, help="Additional args passed through to snakemake")
    return parser.parse_args()


def check_dependencies(include_rscape: bool = False) -> int:
    dependencies = list(DEFAULT_RUNTIME_DEPENDENCIES)
    if include_rscape:
        dependencies.append("R-scape")
    missing: list[str] = []
    for dep in dependencies:
        if shutil.which(dep) is None:
            missing.append(dep)

    if missing:
        print("Missing runtime dependencies:", file=sys.stderr)
        for dep in missing:
            print(f"  - {dep}", file=sys.stderr)
        return 2

    print("All external runtime dependencies are available.")
    return 0


def extract_configfile_arg(snakemake_args: list[str]) -> str | None:
    passthrough = list(snakemake_args)
    if passthrough and passthrough[0] == "--":
        passthrough = passthrough[1:]
    for idx, arg in enumerate(passthrough):
        if arg == "--configfile" and idx + 1 < len(passthrough):
            return passthrough[idx + 1]
        if arg == "--configfiles" and idx + 1 < len(passthrough):
            return passthrough[idx + 1]
    return None


def load_configured_tools(configfile: str) -> dict[str, str]:
    config_path = Path(configfile)
    if not config_path.is_file():
        return {}
    payload = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    tools = payload.get("tools", {}) or {}
    return {str(key): str(value) for key, value in tools.items()}


def show_tool_paths(configfile: str, include_rscape: bool) -> int:
    configured_tools = load_configured_tools(configfile)
    tool_names = [
        "rnalalifold",
        "remove_gaponly",
        "strip_aln",
        "legacy_postprocess",
        "eslreformat",
        "rnaz",
        "alifoldz",
        "rnaalifold",
        "ps2eps",
        "epstopdf",
        "refold",
        "rnafold",
        "magick",
        "cmbuild",
        "cmcalibrate",
    ]
    if include_rscape:
        tool_names.append("rscape")

    default_commands = {
        "rnalalifold": "RNALalifold",
        "remove_gaponly": "python3 -m rnaconsnake.tools.remove_gaponly",
        "strip_aln": "python3 -m rnaconsnake.tools.strip_aln",
        "legacy_postprocess": "python3 -m rnaconsnake.tools.legacy_postprocess",
        "eslreformat": "esl-reformat",
        "rnaz": "RNAz",
        "alifoldz": "alifoldz.pl",
        "rnaalifold": "RNAalifold",
        "ps2eps": "ps2eps",
        "epstopdf": "epstopdf",
        "refold": "refold.pl",
        "rnafold": "RNAfold",
        "magick": "magick",
        "cmbuild": "cmbuild",
        "cmcalibrate": "cmcalibrate",
        "rscape": "R-scape",
    }

    print(f"Tool resolution using config: {Path(configfile).resolve()}")
    for name in tool_names:
        command = configured_tools.get(name, default_commands[name])
        tokens = shlex.split(command)
        executable = tokens[0] if tokens else ""
        resolved = shutil.which(executable) if executable else None
        resolved_text = resolved or "MISSING"
        print(f"{name}:")
        print(f"  command: {command}")
        print(f"  executable: {resolved_text}")
    return 0


def resolve_workflow_paths(snakefile_override: str | None) -> tuple[str, str]:
    if snakefile_override:
        snakefile = Path(snakefile_override).resolve()
        configfile = snakefile.with_name("config.yaml")
        return str(snakefile), str(configfile)

    packaged_root = files("rnaconsnake").joinpath("workflow")
    packaged_snakefile = packaged_root.joinpath("snakefile")
    packaged_config = packaged_root.joinpath("config.yaml")
    if packaged_snakefile.is_file() and packaged_config.is_file():
        return str(packaged_snakefile), str(packaged_config)

    repo_root = Path(__file__).resolve().parents[2]
    repo_snakefile = repo_root / "snakefile"
    repo_config = repo_root / "config.yaml"
    if repo_snakefile.is_file() and repo_config.is_file():
        return str(repo_snakefile), str(repo_config)

    raise FileNotFoundError("Could not locate packaged or repository workflow files.")


JOBSTAT_LINE = re.compile(r"^([A-Za-z0-9_][A-Za-z0-9_.-]*)\s+(\d+)$")
FINISHED_LINE = re.compile(r"Finished jobid:\s+\d+\s+\(Rule:\s+([^)]+)\)")


def should_use_progress(args: argparse.Namespace) -> bool:
    if args.no_progress:
        return False
    if args.progress:
        return tqdm is not None
    return tqdm is not None and sys.stderr.isatty()


def parse_job_stats_lines(lines: list[str]) -> tuple[int, dict[str, int]]:
    total = 0
    counts: dict[str, int] = {}
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        match = JOBSTAT_LINE.match(stripped)
        if not match:
            continue
        name, count_text = match.groups()
        count = int(count_text)
        if name == "total":
            total = count
        elif name not in {"all", "clean"}:
            counts[name] = count
    return total, counts


def run_with_progress(cmd: list[str], env: dict[str, str]) -> int:
    process = subprocess.Popen(
        cmd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    overall = None
    bars: dict[str, tqdm] = {}
    current_stats: list[str] | None = None
    try:
        if process.stdout is None:
            raise RuntimeError("subprocess stdout pipe was not established")
        for raw_line in process.stdout:
            line = raw_line.rstrip("\n")

            if current_stats is not None:
                if not line.strip():
                    total, counts = parse_job_stats_lines(current_stats)
                    if total and overall is None:
                        overall = tqdm(total=total, desc="Overall", position=0, leave=True, dynamic_ncols=True)
                        for position, (rule, count) in enumerate(sorted(counts.items()), start=1):
                            bars[rule] = tqdm(
                                total=count,
                                desc=rule,
                                position=position,
                                leave=True,
                                dynamic_ncols=True,
                            )
                    current_stats = None
                    continue
                current_stats.append(line)
                continue

            if line.strip() == "Job stats:":
                current_stats = []
                continue

            finished = FINISHED_LINE.search(line)
            if finished:
                rule = finished.group(1)
                if overall is not None:
                    overall.update(1)
                bar = bars.get(rule)
                if bar is not None:
                    bar.update(1)
                continue

            if overall is not None:
                tqdm.write(line)
            else:
                print(line)
    finally:
        returncode = process.wait()
        if current_stats is not None:
            total, counts = parse_job_stats_lines(current_stats)
            if total and overall is None:
                overall = tqdm(total=total, desc="Overall", position=0, leave=True, dynamic_ncols=True)
                for position, (rule, count) in enumerate(sorted(counts.items()), start=1):
                    bars[rule] = tqdm(total=count, desc=rule, position=position, leave=True, dynamic_ncols=True)
        if overall is not None:
            overall.close()
        for bar in bars.values():
            bar.close()
    return returncode


def main() -> int:
    args = parse_args()
    if args.version:
        print(f"rnaconsnake {__version__}")
        return 0
    snakefile, default_configfile = resolve_workflow_paths(args.snakefile)
    configfile = extract_configfile_arg(args.snakemake_args) or default_configfile
    if args.show_tool_paths:
        return show_tool_paths(configfile=configfile, include_rscape=args.rscape)
    if args.check_deps:
        return check_dependencies(include_rscape=args.rscape)
    if not args.input_alignment:
        print("Missing required --input-alignment /path/to/input_alignment.{stk,aln}", file=sys.stderr)
        return 2
    dep_status = check_dependencies(include_rscape=args.rscape)
    if dep_status != 0:
        return dep_status

    env = os.environ.copy()
    env.setdefault("PYTHONWARNINGS", "ignore:invalid escape sequence:SyntaxWarning")
    env.setdefault("XDG_CACHE_HOME", tempfile.gettempdir())
    env.setdefault("TMPDIR", tempfile.gettempdir())

    cmd = [
        "snakemake",
        "--snakefile",
        snakefile,
        "--cores",
        str(args.cores),
        "--rerun-incomplete",
        "-p",
        "--latency-wait",
        str(args.latency_wait),
    ]
    if args.output_dir:
        cmd.extend(["--directory", str(Path(args.output_dir).resolve())])
    if args.conservative:
        cmd.extend(["--rerun-triggers", "mtime"])
    if args.snakemake_args:
        passthrough = args.snakemake_args
        if passthrough and passthrough[0] == "--":
            passthrough = passthrough[1:]
        cmd.extend(passthrough)
    config_overrides = [f"input_alignment={Path(args.input_alignment).resolve()}"]
    if args.maxbpspan:
        config_overrides.append(f"maxbpspan=[{','.join(str(value) for value in args.maxbpspan)}]")
    if args.rscape:
        config_overrides.append("do_rscape=True")
    if args.rnaz_shuffle:
        config_overrides.append("rnaz_no_shuffle=False")
    cmd.extend(["--config", *config_overrides])
    if "--configfile" not in cmd and "--configfiles" not in cmd:
        cmd.extend(["--configfile", default_configfile])

    if should_use_progress(args):
        returncode = run_with_progress(cmd, env)
    else:
        returncode = subprocess.call(cmd, env=env)

    if returncode != 0:
        return returncode

    if args.export_bundle:
        run_dir = Path(args.output_dir).resolve() if args.output_dir else Path.cwd().resolve()
        export_args = argparse.Namespace(
            run_dir=str(run_dir),
            output_dir=str(Path(args.export_bundle).resolve()),
            dataset_id=args.export_dataset_id,
            dataset_label=args.export_dataset_label,
            feature_id=args.export_feature_id,
            feature_label=args.export_feature_label,
            feature_type=args.export_feature_type,
            input_alignment=str(Path(args.input_alignment).resolve()),
            source_label=args.export_source_label,
            description=args.export_description,
            overwrite=True,
        )
        build_export(export_args)
        print(f"RNAConSnake export bundle written to: {Path(args.export_bundle).resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
