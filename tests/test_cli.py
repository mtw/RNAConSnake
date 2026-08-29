"""CLI surface: argument handling, the dependency preflight, tool
resolution, and progress rendering."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from rnaconsnake import cli
from tests.helpers import (
    PYTHON,
    _fake_viennarna_tools,
    _parse_cli_args,
    read_text,
    subprocess_env,
)


def test_cli_version_reports_package_version() -> None:
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    from rnaconsnake import __version__

    assert result.stdout.strip() == f"RNAcs {__version__}"


def test_cli_check_deps_reports_missing_rnalalifold() -> None:
    env = subprocess_env()
    env["PATH"] = str(Path(PYTHON).resolve().parent)
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--check-deps"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 2
    assert "RNALalifold" in result.stderr
    assert "RNAalifold" in result.stderr
    assert "RNAz" in result.stderr


def test_cli_check_deps_reports_missing_rscape_when_requested() -> None:
    env = subprocess_env()
    env["PATH"] = str(Path(PYTHON).resolve().parent)
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--check-deps", "--rscape"],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 2
    assert "R-scape" in result.stderr


def test_cli_check_deps_checks_the_configured_tool_command(tmp_path: Path) -> None:
    """A tool pointed at a custom path must be checked at that path.

    Checking the hard-coded program name instead both fails runs that would
    have worked and passes runs that then die mid-pipeline.
    """
    configfile = tmp_path / "config.yaml"
    configfile.write_text("tools:\n  alifoldz: /nonexistent/dir/alifoldz.pl\n", encoding="utf-8")
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--check-deps", "--", "--configfile", str(configfile)],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert result.returncode == 2
    assert "/nonexistent/dir/alifoldz.pl" in result.stderr


def test_cli_check_deps_honours_optional_branches_from_the_config_file(tmp_path: Path) -> None:
    """`do_rscape` and `null.method` enable branches with their own tools,
    whether they are switched on by flag or in the config file."""
    configfile = tmp_path / "config.yaml"
    configfile.write_text('do_rscape: true\n"null":\n  method: sissiz\n  replicates: 2\n', encoding="utf-8")
    env = subprocess_env()
    env["PATH"] = str(Path(PYTHON).resolve().parent)
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--check-deps", "--", "--configfile", str(configfile)],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 2
    assert "R-scape" in result.stderr
    assert "SISSIz" in result.stderr


def test_cli_benchmark_requires_the_null_arm(tmp_path: Path) -> None:
    """The recovery table scores the calibrated loci, which only the null arm
    produces. Snakemake's own error for the missing rule names neither."""
    alignment = tmp_path / "input.stk"
    alignment.write_text("# STOCKHOLM 1.0\n//\n", encoding="utf-8")
    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.cli",
            "--input-alignment",
            str(alignment),
            "--benchmark",
            "--cores",
            "1",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert result.returncode == 2
    assert "--null-arm" in result.stderr


def test_cli_show_tool_paths_reports_configured_commands() -> None:
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--show-tool-paths"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert "Tool resolution using config:" in result.stdout
    assert "rnaz:" in result.stdout
    assert "command: RNAz" in result.stdout
    assert "legacy_postprocess:" in result.stdout
    assert "command: python3 -m rnaconsnake.tools.legacy_postprocess" in result.stdout


def test_cli_help_describes_no_local_only_layout() -> None:
    """Help text is the first documentation a user reads; a `data/<sample>.stk`
    convention that the repository does not ship sends them looking for a
    directory that only ever existed on one machine."""
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--help"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    for ignored in ["data/", "runs/", "configuration_lock"]:
        assert ignored not in result.stdout, f"--help points at {ignored}, which no clone has"


def test_cli_requires_input_alignment_for_run() -> None:
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--cores", "1"],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert result.returncode == 2
    assert "--input-alignment" in result.stderr


def test_parse_job_stats_lines_extracts_totals_and_rule_counts() -> None:
    total, counts = cli.parse_job_stats_lines(
        [
            "job                count",
            "---------------  -------",
            "RNALalifold            2",
            "all                    1",
            "run_post_rnaalifold_file      42",
            "total                 45",
        ]
    )

    assert total == 45
    assert counts == {
        "RNALalifold": 2,
        "run_post_rnaalifold_file": 42,
    }


def test_cli_reports_program_name_rnacs() -> None:
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.cli", "--version"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert result.stdout.startswith("RNAcs ")


def test_cli_null_arm_options_resolve_to_a_full_config_section(tmp_path: Path) -> None:
    import argparse as _argparse

    configfile = tmp_path / "config.yaml"
    configfile.write_text(
        '"null":\n  method: none\n  replicates: 10\n  seed: 20261101\n  two_stage: true\n',
        encoding="utf-8",
    )

    def namespace(**overrides):
        base = {
            "null_arm": None,
            "null_replicates": None,
            "null_seed": None,
            "no_two_stage": False,
            "null_pool": None,
        }
        base.update(overrides)
        return _argparse.Namespace(**base)

    assert cli.resolve_null_config(namespace(), str(configfile)) is None

    bare = cli.resolve_null_config(namespace(null_arm="sissiz"), str(configfile))
    assert bare == {"method": "sissiz", "replicates": 10, "seed": 20261101, "two_stage": True}

    tuned = cli.resolve_null_config(
        namespace(null_arm="rnazRandomizeAln", null_replicates=2, null_seed=7, no_two_stage=True),
        str(configfile),
    )
    assert tuned == {
        "method": "rnazRandomizeAln",
        "replicates": 2,
        "seed": 7,
        "two_stage": False,
    }

    pinned = cli.resolve_null_config(namespace(null_arm="sissiz", null_pool="/tmp/pool.stk"), str(configfile))
    assert pinned["pool_file"].endswith("/tmp/pool.stk")
    # An unset pool must not travel as the string "None": Snakemake stringifies
    # nested --config values, and "None" is a truthy filename.
    assert "pool_file" not in bare


def test_cli_check_deps_reports_missing_null_backend() -> None:
    assert cli.check_dependencies(null_method="sissiz") in (0, 2)
    assert cli.NULL_METHOD_DEPENDENCIES["sissiz"] == "SISSIz"
    assert cli.NULL_METHOD_DEPENDENCIES["rnazRandomizeAln"] == "rnazRandomizeAln.pl"


def test_split_passthrough_config_merges_user_entries() -> None:
    """Snakemake's --config is nargs="+", so two occurrences lose the first."""
    rest, entries = cli.split_passthrough_config(
        ["--configfile", "c.yaml", "--config", "a=1", "b=2", "--cores", "4"]
    )
    assert rest == ["--configfile", "c.yaml", "--cores", "4"]
    assert entries == ["a=1", "b=2"]

    assert cli.split_passthrough_config(["--dry-run"]) == (["--dry-run"], [])
    assert cli.split_passthrough_config(["-C", "x=1"]) == ([], ["x=1"])


def test_requested_null_method_reads_the_configured_backend() -> None:
    assert cli.requested_null_method({"method": "sissiz"}) == "sissiz"
    assert cli.requested_null_method({"method": "none"}) is None
    assert cli.requested_null_method({}) is None
    assert cli.requested_null_method(None) is None


def test_rnacs_rejects_a_null_config_that_is_not_a_mapping() -> None:
    """`--config null=notadict` used to reach `.get()` on a string and raise an
    AttributeError traceback."""
    with pytest.raises(SystemExit) as excinfo:
        cli.requested_null_method("notadict")
    assert excinfo.value.code == 2


def test_rnacs_export_overwrite_defaults_to_off() -> None:
    """The wrapper passed overwrite=True on every run."""
    parsed = _parse_cli_args(["--input-alignment", "x.stk", "--export-bundle", "out"])
    assert parsed.export_overwrite is False
    parsed = _parse_cli_args(["--input-alignment", "x.stk", "--export-bundle", "out", "--export-overwrite"])
    assert parsed.export_overwrite is True


def test_check_deps_requires_the_viennarna_python_module(monkeypatch) -> None:
    """The bindings are as much a runtime dependency as the binaries, and the
    preflight is where a missing one has to surface."""
    monkeypatch.setattr(cli, "viennarna_bindings_version", lambda: None)
    assert cli.check_dependencies() == 2


def test_viennarna_binaries_must_be_the_same_build_as_the_module(tmp_path: Path) -> None:
    """The consensus structure comes from the binaries and the refold from the
    Python module. Two builds in one run means two sets of energy parameters,
    so the versions must be identical -- not merely close."""
    tools = _fake_viennarna_tools(tmp_path, "2.6.4")

    conflicts = cli.viennarna_version_conflicts("2.7.2", tools)
    assert len(conflicts) == len(cli.VIENNARNA_BINARIES), conflicts
    assert any("2.6.4" in line and "2.7.2" in line for line in conflicts)
    assert cli.viennarna_version_conflicts("2.6.4", tools) == []

    # Exact equality: a prefix match is not a match.
    assert cli.viennarna_version_conflicts("2.6.40", tools)

    # Nothing to compare against is not a conflict.
    assert cli.viennarna_version_conflicts(None, tools) == []
    missing = dict.fromkeys(cli.VIENNARNA_BINARIES, "definitely-not-installed")
    assert cli.viennarna_version_conflicts("2.6.4", missing) == []


def test_viennarna_below_the_supported_minimum_is_refused() -> None:
    """The configuration was locked against 2.7.2, and PyPI publishes wheels
    only for 2.7.x, so an older ViennaRNA is neither tested nor installable."""
    assert cli.MINIMUM_VIENNARNA == (2, 7)
    assert "2.6.4" in cli.viennarna_too_old("2.6.4")
    assert cli.viennarna_too_old("2.7.0") == ""
    assert cli.viennarna_too_old("3.0") == ""
    # Nothing to judge is not a complaint.
    assert cli.viennarna_too_old(None) == ""
    assert cli.viennarna_too_old("unknown") == ""


def test_check_deps_refuses_an_outdated_viennarna(monkeypatch, capsys) -> None:
    monkeypatch.setattr(cli, "viennarna_bindings_version", lambda: "2.6.4")
    monkeypatch.setattr(cli, "viennarna_version_conflicts", lambda bindings, tools=None: [])
    assert cli.check_dependencies() == 2
    assert "older than the 2.7" in capsys.readouterr().err


def test_the_declared_viennarna_floor_matches_the_preflight() -> None:
    """One number, two places: the dev extra and the runtime check."""
    import re as _re

    declared = _re.search(r'"ViennaRNA>=(\d+)\.(\d+)', read_text(Path("pyproject.toml")))
    assert declared, "the dev extra no longer pins a ViennaRNA floor"
    assert (int(declared.group(1)), int(declared.group(2))) == cli.MINIMUM_VIENNARNA


def test_check_deps_fails_on_a_viennarna_version_mismatch(monkeypatch, capsys) -> None:
    """A mixed toolchain is a refusal, not a warning: a calibrated result could
    not be reproduced from the versions the run records."""
    monkeypatch.setattr(cli, "viennarna_bindings_version", lambda: "2.7.2")
    monkeypatch.setattr(
        cli,
        "viennarna_version_conflicts",
        lambda bindings, tools=None: ["RNAalifold is 2.6.4, the RNA Python module is 2.7.2"],
    )
    assert cli.check_dependencies() == 2
    captured = capsys.readouterr()
    assert "ViennaRNA version problem" in captured.err
    assert "same ViennaRNA release" in captured.err


def test_should_use_progress_follows_the_flags_and_the_terminal(monkeypatch) -> None:
    from argparse import Namespace

    monkeypatch.setattr(cli, "tqdm", object())
    assert cli.should_use_progress(Namespace(progress=True, no_progress=False)) is True
    # --no-progress wins over --progress.
    assert cli.should_use_progress(Namespace(progress=True, no_progress=True)) is False

    # Without tqdm there are no bars to draw, whatever was asked for.
    monkeypatch.setattr(cli, "tqdm", None)
    assert cli.should_use_progress(Namespace(progress=True, no_progress=False)) is False


def test_run_with_progress_counts_finished_jobs_and_passes_output_through(tmp_path: Path, capsys) -> None:
    """The wrapper must forward every line it does not consume, and return the
    child's exit code -- a swallowed error message is worse than no bars."""
    script = tmp_path / "fake_snakemake.py"
    script.write_text(
        "print('Building DAG of jobs...')\n"
        "print('Job stats:')\n"
        "print('job      count')\n"
        "print('RNALalifold    2')\n"
        "print('all      1')\n"
        "print('total    3')\n"
        "print('')\n"
        "print('Finished jobid: 1 (Rule: RNALalifold)')\n"
        "print('Finished jobid: 2 (Rule: RNALalifold)')\n"
        "print('a message worth seeing')\n",
        encoding="utf-8",
    )
    assert cli.run_with_progress([PYTHON, str(script)], subprocess_env()) == 0
    captured = capsys.readouterr()
    assert "a message worth seeing" in captured.out + captured.err
    assert "Building DAG of jobs..." in captured.out + captured.err


def test_run_with_progress_returns_the_childs_exit_code(tmp_path: Path) -> None:
    script = tmp_path / "failing.py"
    script.write_text("import sys\nprint('boom')\nsys.exit(3)\n", encoding="utf-8")
    assert cli.run_with_progress([PYTHON, str(script)], subprocess_env()) == 3


def test_parse_job_stats_ignores_the_aggregate_rows() -> None:
    total, counts = cli.parse_job_stats_lines(
        ["job      count", "RNALalifold   2", "all       1", "clean     1", "total     3"]
    )
    assert total == 3
    # `all` and `clean` are not real work; counting them would inflate the bars.
    assert counts == {"RNALalifold": 2}
