from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from importlib.resources import files
from pathlib import Path

from rnaconsnake import __version__

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
        "--output-dir",
        help="Working directory where all workflow outputs will be written.",
    )
    parser.add_argument(
        "--input-alignment",
        help="Run the workflow on an arbitrary input Stockholm alignment file instead of data/<sample>.stk.",
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
        assert process.stdout is not None
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
    if args.check_deps:
        return check_dependencies(include_rscape=args.rscape)
    if not args.input_alignment:
        print("Missing required --input-alignment /path/to/input.stk", file=sys.stderr)
        return 2
    dep_status = check_dependencies(include_rscape=args.rscape)
    if dep_status != 0:
        return dep_status

    snakefile, default_configfile = resolve_workflow_paths(args.snakefile)
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
    cmd.extend(["--config", *config_overrides])
    if "--configfile" not in cmd and "--configfiles" not in cmd:
        cmd.extend(["--configfile", default_configfile])

    if should_use_progress(args):
        return run_with_progress(cmd, env)
    return subprocess.call(cmd, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
