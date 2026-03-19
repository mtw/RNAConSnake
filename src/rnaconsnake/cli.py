from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from importlib.resources import files
from pathlib import Path

from rnaconsnake import __version__


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the packaged RNAConSnake Snakemake workflow.")
    parser.add_argument("--cores", default="all", help="Snakemake cores value, e.g. 4 or all")
    parser.add_argument("--latency-wait", default="20", help="Snakemake latency wait")
    parser.add_argument("--snakefile", help="Override packaged snakefile path")
    parser.add_argument("--check-deps", action="store_true", help="Check external runtime dependencies and exit")
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    parser.add_argument("snakemake_args", nargs=argparse.REMAINDER, help="Additional args passed through to snakemake")
    return parser.parse_args()


def check_dependencies() -> int:
    missing: list[str] = []
    if shutil.which("RNALalifold") is None:
        missing.append("RNALalifold")
    if shutil.which("snakemake") is None:
        missing.append("snakemake")

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


def main() -> int:
    args = parse_args()
    if args.version:
        print(f"rnaconsnake {__version__}")
        return 0
    if args.check_deps:
        return check_dependencies()
    dep_status = check_dependencies()
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
    if args.snakemake_args:
        passthrough = args.snakemake_args
        if passthrough and passthrough[0] == "--":
            passthrough = passthrough[1:]
        cmd.extend(passthrough)
    if "--configfile" not in cmd and "--configfiles" not in cmd:
        cmd.extend(["--configfile", default_configfile])

    return subprocess.call(cmd, env=env)


if __name__ == "__main__":
    raise SystemExit(main())
