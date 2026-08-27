from __future__ import annotations

import argparse
import ast
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
from rnaconsnake.workflow_helpers import NullSettings

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback only
    tqdm = None

# Every external command the workflow can invoke, and the command used when
# ``tools:`` in the config file does not override it. The preflight check and
# ``--show-tool-paths`` both resolve through this table, so a configured
# ``tools.alifoldz: /custom/path/alifoldz.pl`` is what actually gets checked.
DEFAULT_TOOL_COMMANDS = {
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
    "refold": "python3 -m rnaconsnake.tools.refold",
    "magick": "magick",
    "cmbuild": "cmbuild",
    "cmcalibrate": "cmcalibrate",
    "rscape": "R-scape",
    "sissiz": "SISSIz",
    "rnaz_randomize_aln": "rnazRandomizeAln.pl",
    "null_model": "python3 -m rnaconsnake.tools.null_model",
    "calibration": "python3 -m rnaconsnake.tools.calibration",
    "benchmark": "python3 -m rnaconsnake.tools.benchmark",
    "dereplicate": "python3 -m rnaconsnake.tools.dereplicate",
    "versions": "python3 -m rnaconsnake.tools.versions",
    "alignment_report": "python3 -m rnaconsnake.tools.alignment_report",
    "benchmark_scaffold": "python3 -m rnaconsnake.tools.benchmark_scaffold",
    "threshold_sweep": "python3 -m rnaconsnake.tools.threshold_sweep",
    "alifoldmaxcovar": "python3 -m rnaconsnake.tools.alifold_maxcovar",
}

# Tool keys every run needs, regardless of which optional branches are on.
DEFAULT_RUNTIME_TOOLS = [
    "rnalalifold",
    "eslreformat",
    "rnaz",
    "alifoldz",
    "rnaalifold",
    "ps2eps",
    "epstopdf",
    "magick",
]

# Not a configurable tool: RNAcs invokes it directly.
DEFAULT_RUNTIME_DEPENDENCIES = ["snakemake"]


PROGRAM_NAME = "RNAcs"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog=PROGRAM_NAME,
        description="Run the packaged RNAConSnake Snakemake workflow.",
    )
    parser.add_argument("--cores", default="all", help="Snakemake cores value, e.g. 4 or all")
    parser.add_argument("--latency-wait", default="20", help="Snakemake latency wait")
    parser.add_argument("--snakefile", help="Override packaged snakefile path")
    parser.add_argument(
        "--maxbpspan",
        type=int,
        action="append",
        help=(
            "Override maxbpspan values from config.yaml. Repeat the option to provide multiple window sizes."
        ),
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
        help="Input alignment to screen, in Stockholm (.stk) or Clustal (.aln) format. Required.",
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
        "--null-arm",
        nargs="?",
        const="sissiz",
        choices=["sissiz", "rnazRandomizeAln", "none"],
        help=(
            "Enable the null-model calibration arm. Bare --null-arm uses SISSIz, the method "
            "the manuscript reports. Every pipeline output then moves under arms/<arm>/ and "
            "results/calibration/ is produced."
        ),
    )
    parser.add_argument(
        "--null-replicates",
        type=int,
        help="Number of null replicates (10 = sanity check, 100 = usable q-values).",
    )
    parser.add_argument(
        "--null-seed",
        type=int,
        help="Base seed for null replicate generation.",
    )
    parser.add_argument(
        "--null-pool",
        help=(
            "Reuse a previously generated null_pool/pool.stk instead of simulating new "
            "replicates. SISSIz cannot be seeded, so pinning the pool is how a "
            "calibration is reproduced."
        ),
    )
    parser.add_argument(
        "--no-two-stage",
        action="store_true",
        help=(
            "Run AlifoldZ on every candidate instead of only on stage-one survivors. "
            "Much slower, but the resulting FDR is unconditional."
        ),
    )
    parser.add_argument(
        "--dereplicate",
        choices=["containment", "substructure", "overlap", "none"],
        help=(
            "How to collapse overlapping RNALalifold windows into one candidate per locus. "
            "Default: containment (a window nested inside another is a fragment of it)."
        ),
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Also build the flavivirus positive-control recovery table (needs a curated truth file).",
    )
    parser.add_argument(
        "--emit-versions",
        action="store_true",
        help="Write results/versions.yaml recording the toolchain used for this run.",
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
    parser.add_argument(
        "--check-deps", action="store_true", help="Check external runtime dependencies and exit"
    )
    parser.add_argument("--version", action="store_true", help="Show version and exit")
    parser.add_argument(
        "snakemake_args", nargs=argparse.REMAINDER, help="Additional args passed through to snakemake"
    )
    return parser.parse_args()


# ViennaRNA binaries: every one of these comes from the same install as the
# Python module, and must report the same version as it.
VIENNARNA_BINARIES = ("rnalalifold", "rnaalifold")

# The oldest ViennaRNA this runs on. The shipped configuration was locked
# against 2.7.2, and PyPI publishes wheels only for 2.7.x, so an older bound
# would be one nothing can install anyway.
MINIMUM_VIENNARNA = (2, 7)

VERSION_NUMBER = re.compile(r"\d+\.\d+(?:\.\d+)*")


NULL_METHOD_TOOLS = {
    "sissiz": "sissiz",
    "rnazRandomizeAln": "rnaz_randomize_aln",
}

# The default executable behind each null-model backend, for messages that name
# the program rather than the tool key.
NULL_METHOD_DEPENDENCIES = {method: DEFAULT_TOOL_COMMANDS[name] for method, name in NULL_METHOD_TOOLS.items()}


def tool_command(name: str, configured_tools: dict[str, str] | None = None) -> str:
    """The command configured for a tool, falling back to its default."""
    configured = (configured_tools or {}).get(name)
    return str(configured) if configured else DEFAULT_TOOL_COMMANDS[name]


def required_tool_names(include_rscape: bool = False, null_method: str | None = None) -> list[str]:
    tools = list(DEFAULT_RUNTIME_TOOLS)
    if include_rscape:
        tools.append("rscape")
    if null_method in NULL_METHOD_TOOLS:
        tools.append(NULL_METHOD_TOOLS[null_method])
    return tools


def viennarna_bindings_version() -> str | None:
    """``RNA.__version__``, or ``None`` when the module is not importable.

    Refolding runs through the bindings rather than ``refold.pl``, so they are
    as much a runtime dependency as the binaries are.
    """
    try:
        import RNA
    except ImportError:
        return None
    return str(getattr(RNA, "__version__", "unknown"))


def check_dependencies(
    include_rscape: bool = False,
    null_method: str | None = None,
    configured_tools: dict[str, str] | None = None,
) -> int:
    """Check that every command this run will invoke actually resolves.

    The commands come from the run's own ``tools:`` configuration, so a
    tool pointed at a custom path passes, and a tool pointed at a path that
    does not exist fails here rather than mid-run.
    """
    missing: list[str] = []
    for dep in DEFAULT_RUNTIME_DEPENDENCIES:
        if shutil.which(dep) is None:
            missing.append(dep)
    bindings = viennarna_bindings_version()
    if bindings is None:
        missing.append("the ViennaRNA Python module (conda: viennarna, pip: ViennaRNA)")
    for name in required_tool_names(include_rscape=include_rscape, null_method=null_method):
        command = tool_command(name, configured_tools)
        tokens = shlex.split(command)
        executable = tokens[0] if tokens else ""
        if not executable or shutil.which(executable) is None:
            missing.append(command if command == executable else f"{executable} (from {name}: {command})")

    conflicts = viennarna_version_conflicts(bindings, configured_tools)
    outdated = viennarna_too_old(bindings)
    if outdated:
        conflicts.insert(0, outdated)

    if missing or conflicts:
        if missing:
            print("Missing runtime dependencies:", file=sys.stderr)
            for dep in missing:
                print(f"  - {dep}", file=sys.stderr)
        if conflicts:
            print("ViennaRNA version problem:", file=sys.stderr)
            for conflict in conflicts:
                print(f"  - {conflict}", file=sys.stderr)
            print(
                "  The consensus structure comes from the binaries and the refold from "
                "the module.\n"
                "  Install both from the same ViennaRNA release; one run must not mix "
                "two builds'\n  energy parameters.",
                file=sys.stderr,
            )
        return 2

    print(f"All external runtime dependencies are available (ViennaRNA {parse_version(bindings)}).")
    return 0


def viennarna_version_conflicts(
    bindings: str | None, configured_tools: dict[str, str] | None = None
) -> list[str]:
    """ViennaRNA binaries whose version differs from the Python module's.

    The consensus structure comes from the binaries and the refold from the
    module. Two ViennaRNA builds in one run means two sets of energy
    parameters, so a calibrated result could not be reproduced from the
    toolchain it records. The versions have to be identical, not merely close.
    """
    wanted = parse_version(bindings)
    if not wanted:
        return []
    conflicts = []
    for name in VIENNARNA_BINARIES:
        command = tool_command(name, configured_tools)
        reported = probe_version(command)
        found = parse_version(reported)
        if found and found != wanted:
            conflicts.append(f"{command} is {found}, the RNA Python module is {wanted}")
    return conflicts


def parse_version(text: str | None) -> str:
    """The first dotted version number in ``text`` (``RNAalifold 2.7.2`` -> ``2.7.2``)."""
    match = VERSION_NUMBER.search(text or "")
    return match.group(0) if match else ""


def version_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split(".")) if version else ()


def viennarna_too_old(bindings: str | None) -> str:
    """Complaint about a ViennaRNA older than RNAcs supports, or ``""``."""
    found = version_tuple(parse_version(bindings))
    if not found or found >= MINIMUM_VIENNARNA:
        return ""
    minimum = ".".join(str(part) for part in MINIMUM_VIENNARNA)
    return f"ViennaRNA {parse_version(bindings)} is older than the {minimum} this needs"


def probe_version(command: str) -> str:
    """First line of ``<command> --version``, empty when it cannot be asked."""
    tokens = shlex.split(command)
    if not tokens or shutil.which(tokens[0]) is None:
        return ""
    try:
        # In a scratch directory: ViennaRNA tools write side outputs into
        # whatever directory they are run from, and this one runs wherever the
        # user happened to invoke RNAcs.
        with tempfile.TemporaryDirectory() as scratch:
            result = subprocess.run(
                [*tokens, "--version"],
                stdin=subprocess.DEVNULL,
                capture_output=True,
                text=True,
                check=False,
                timeout=20,
                cwd=scratch,
            )
    except (OSError, subprocess.SubprocessError):
        return ""
    for line in ((result.stdout or "") + "\n" + (result.stderr or "")).splitlines():
        if line.strip():
            return line.strip()
    return ""


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


def load_config_payload(configfile: str) -> dict:
    config_path = Path(configfile)
    if not config_path.is_file():
        return {}
    return yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}


def load_configured_tools(configfile: str) -> dict[str, str]:
    tools = load_config_payload(configfile).get("tools", {}) or {}
    return {str(key): str(value) for key, value in tools.items()}


def show_tool_paths(configfile: str, include_rscape: bool) -> int:
    configured_tools = load_configured_tools(configfile)
    # Union, so a tool key configured for this run is reported even if RNAcs
    # itself has no default for it.
    tool_names = sorted((set(DEFAULT_TOOL_COMMANDS) | set(configured_tools)) - {"rscape"})
    if include_rscape:
        tool_names.append("rscape")

    print(f"Tool resolution using config: {Path(configfile).resolve()}")
    for name in tool_names:
        command = tool_command(name, configured_tools)
        tokens = shlex.split(command)
        executable = tokens[0] if tokens else ""
        resolved = shutil.which(executable) if executable else None
        resolved_text = resolved or "MISSING"
        print(f"{name}:")
        print(f"  command: {command}")
        print(f"  executable: {resolved_text}")
    return 0


def split_passthrough_config(passthrough: list[str]) -> tuple[list[str], list[str]]:
    """Pull ``--config`` entries out of the passthrough arguments.

    Snakemake's ``--config`` takes ``nargs="+"``, so argparse keeps only the
    last occurrence. Emitting our own ``--config`` alongside the user's would
    silently drop theirs, so both sets are merged into a single option.
    """
    rest: list[str] = []
    entries: list[str] = []
    index = 0
    while index < len(passthrough):
        arg = passthrough[index]
        if arg not in {"--config", "-C"}:
            rest.append(arg)
            index += 1
            continue
        index += 1
        while index < len(passthrough) and not passthrough[index].startswith("-"):
            entries.append(passthrough[index])
            index += 1
    return rest, entries


def load_config_section(configfile: str, key: str) -> dict:
    payload = load_config_payload(configfile)
    # An unquoted ``null:`` key in YAML parses as the null scalar.
    section = payload.get(key)
    if section is None and key == "null":
        section = payload.get(None)
    return dict(section or {})


def passthrough_config_overrides(snakemake_args: list[str]) -> dict:
    """``--config key=value`` entries the user passed through to snakemake.

    The preflight has to see them: ``--config do_rscape=True`` enables a branch
    with its own dependency, exactly as ``--rscape`` does. Values are parsed as
    Python literals where possible, which is how snakemake reads them too.
    """
    passthrough = list(snakemake_args or [])
    if passthrough and passthrough[0] == "--":
        passthrough = passthrough[1:]
    _, entries = split_passthrough_config(passthrough)
    overrides: dict = {}
    for entry in entries:
        if "=" not in entry:
            continue
        key, _, raw = entry.partition("=")
        try:
            overrides[key.strip()] = ast.literal_eval(raw)
        except (SyntaxError, ValueError):
            overrides[key.strip()] = raw
    return overrides


def _truthy(*values) -> bool:
    for value in values:
        if isinstance(value, str):
            if value.strip().lower() in {"true", "yes", "1"}:
                return True
        elif value:
            return True
    return False


def null_arm_enabled(null_settings) -> bool:
    """Whether a run with these ``null:`` settings produces a calibration.

    Decided by the same ``NullSettings`` the workflow itself uses, so the CLI
    cannot disagree with the snakefile about whether a run is calibrated.
    """
    if not isinstance(null_settings, dict):
        return False
    try:
        return NullSettings.from_config({"null": null_settings}).enabled
    except (TypeError, ValueError):
        return False


def resolve_null_config(args: argparse.Namespace, configfile: str) -> dict | None:
    """Merge the null-arm CLI options over the configured defaults.

    Returns ``None`` when no null-arm option was given, so the config file is
    left untouched and a plain run keeps its current behaviour.
    """
    requested = [
        args.null_arm,
        args.null_replicates,
        args.null_seed,
        args.no_two_stage,
        args.null_pool,
    ]
    if not any(value not in (None, False) for value in requested):
        return None

    resolved = {
        "method": "none",
        "replicates": 10,
        "seed": 20261101,
        "two_stage": True,
    }
    resolved.update(load_config_section(configfile, "null"))
    if args.null_arm is not None:
        resolved["method"] = args.null_arm
    elif resolved["method"] == "none":
        # --null-replicates/--null-seed on their own still mean "calibrate".
        resolved["method"] = "sissiz"
    if args.null_replicates is not None:
        resolved["replicates"] = args.null_replicates
    if args.null_seed is not None:
        resolved["seed"] = args.null_seed
    if args.no_two_stage:
        resolved["two_stage"] = False
    if args.null_pool:
        resolved["pool_file"] = str(Path(args.null_pool).resolve())
    # Snakemake stringifies nested --config values, so an empty pool_file would
    # arrive as the literal "None". Omit the key instead.
    if not resolved.get("pool_file"):
        resolved.pop("pool_file", None)
    return resolved


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
                        overall = tqdm(
                            total=total, desc="Overall", position=0, leave=True, dynamic_ncols=True
                        )
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
                    bars[rule] = tqdm(
                        total=count, desc=rule, position=position, leave=True, dynamic_ncols=True
                    )
        if overall is not None:
            overall.close()
        for bar in bars.values():
            bar.close()
    return returncode


def main() -> int:
    args = parse_args()
    if args.version:
        print(f"{PROGRAM_NAME} {__version__}")
        return 0
    snakefile, default_configfile = resolve_workflow_paths(args.snakefile)
    configfile = extract_configfile_arg(args.snakemake_args) or default_configfile
    if args.show_tool_paths:
        return show_tool_paths(configfile=configfile, include_rscape=args.rscape)
    configured_tools = load_configured_tools(configfile)
    overrides = passthrough_config_overrides(args.snakemake_args)
    null_config = resolve_null_config(args, configfile)
    # A run is calibrated when the CLI asks for it, when a --config override
    # asks for it, or when the config file already has the arm switched on.
    effective_null = null_config or overrides.get("null") or load_config_section(configfile, "null")
    null_method = effective_null.get("method") if effective_null else None
    if null_method == "none":
        null_method = None
    include_rscape = args.rscape or _truthy(
        overrides.get("do_rscape"), load_config_payload(configfile).get("do_rscape")
    )
    if args.check_deps:
        return check_dependencies(
            include_rscape=include_rscape,
            null_method=null_method,
            configured_tools=configured_tools,
        )
    if not args.input_alignment:
        print("Missing required --input-alignment /path/to/input_alignment.{stk,aln}", file=sys.stderr)
        return 2
    if args.benchmark and not null_arm_enabled(effective_null):
        # The recovery table scores results/calibration/qvalues.tsv, which only
        # the null arm produces; without it snakemake fails with a missing-rule
        # error that says nothing about the cause.
        print(
            "--benchmark scores the calibrated loci, so it needs the null-model arm. "
            "Add --null-arm sissiz (or set null.method in the config file).",
            file=sys.stderr,
        )
        return 2
    dep_status = check_dependencies(
        include_rscape=include_rscape,
        null_method=null_method,
        configured_tools=configured_tools,
    )
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
    if args.benchmark:
        # Positional targets have to come before any nargs="+" option, ours or
        # the user's, or argparse swallows them as malformed --config entries.
        # "all" is named explicitly: naming any target at all replaces the
        # default one, and --benchmark builds the recovery table *in addition
        # to* the run, not instead of it.
        cmd += ["all", "results/benchmark/flavivirus_recovery.tsv"]
    if args.output_dir:
        cmd.extend(["--directory", str(Path(args.output_dir).resolve())])
    if args.conservative:
        cmd.extend(["--rerun-triggers", "mtime"])
    passthrough_config: list[str] = []
    if args.snakemake_args:
        passthrough = args.snakemake_args
        if passthrough and passthrough[0] == "--":
            passthrough = passthrough[1:]
        passthrough, passthrough_config = split_passthrough_config(passthrough)
        cmd.extend(passthrough)
    # The user's entries come first so that explicit RNAcs flags win.
    config_overrides = [*passthrough_config, f"input_alignment={Path(args.input_alignment).resolve()}"]
    if args.maxbpspan:
        config_overrides.append(f"maxbpspan=[{','.join(str(value) for value in args.maxbpspan)}]")
    if args.rscape:
        config_overrides.append("do_rscape=True")
    if args.rnaz_shuffle:
        config_overrides.append("rnaz_no_shuffle=False")
    if null_config is not None:
        # Snakemake replaces a --config key wholesale rather than merging, so
        # the whole section is passed at once.
        config_overrides.append("null=" + repr(null_config))
    if args.dereplicate:
        config_overrides.append("dereplicate=" + repr({"method": args.dereplicate}))
    if args.emit_versions:
        config_overrides.append("emit_versions=True")
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
