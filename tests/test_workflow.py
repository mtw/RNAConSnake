"""Workflow topology, canonical paths, the Snakemake smoke tests,
and run-to-run consistency."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

from rnaconsnake import cli
from rnaconsnake.workflow_helpers import CandidatePaths, WorkflowSettings, initial_alignment_format_code
from tests.helpers import (
    _REFOLD_RECORD,
    PYTHON,
    ROOT,
    TOY_ALIGNMENT,
    _run_verify,
    _write_run_tree,
    read_text,
    subprocess_env,
    write_fake_postprocess_tools,
    write_fake_randomize_aln,
    write_fake_rnalalifold,
)


def test_candidate_paths_exposes_canonical_workflow_locations() -> None:
    paths = CandidatePaths(wlen=150, file="RC_150_0001_aln_1_12")

    assert paths.split == "Lalifold/len_150/split/RC_150_0001_aln_1_12.stk"
    assert paths.remgap == "generated_files/remgap/len_150/RC_150_0001_aln_1_12_remgap.stk"
    assert paths.strip == "generated_files/strip/len_150/RC_150_0001_aln_1_12_stripped.stk"
    assert paths.stk == "generated_files/stk/len_150/RC_150_0001_aln_1_12.stk"
    assert paths.aln == "generated_files/aln/len_150/RC_150_0001_aln_1_12.aln"
    assert paths.rnaz_json == "generated_files/rnaz/len_150/RC_150_0001_aln_1_12.rnaz.json"
    assert paths.alifoldz_json == "generated_files/alifoldz/len_150/RC_150_0001_aln_1_12.alifoldz.json"
    assert paths.rnaalifold_stk == (
        "generated_files/rnaalifold/len_150/RC_150_0001_aln_1_12/RC_150_0001_aln_1_12.RNAalifold_results.stk"
    )
    assert paths.consensus_json == ("generated_files/consensus/len_150/RC_150_0001_aln_1_12.consensus.json")
    assert paths.maxcovar_json == "generated_files/maxcovar/len_150/RC_150_0001_aln_1_12.maxcovar.json"
    assert paths.rscape_json == "generated_files/rscape/len_150/RC_150_0001_aln_1_12.rscape.json"
    assert paths.summary_json == "generated_files/summary/len_150/RC_150_0001_aln_1_12.summary.json"
    assert paths.png_aln == "generated_files/png/len_150/RC_150_0001_aln_1_12_aln.png"
    assert paths.cm_status_json == "generated_files/cm/len_150/RC_150_0001_aln_1_12.cm.status.json"


def test_initial_alignment_format_code_detects_stockholm_and_clustal() -> None:
    assert initial_alignment_format_code("example.stk") == "S"
    assert initial_alignment_format_code("example.STK") == "S"
    assert initial_alignment_format_code("example.aln") == "C"
    assert initial_alignment_format_code("example.ALN") == "C"
    assert initial_alignment_format_code("example.anything_else") == "S"


def test_initial_alignment_format_code_handles_suffixless_input() -> None:
    assert initial_alignment_format_code("example") == "S"
    assert initial_alignment_format_code("example.") == "S"


def test_initial_alignment_format_code_missing_input() -> None:
    with pytest.raises(ValueError, match=r"input_alignment") as exc:
        initial_alignment_format_code(None)
    assert "input_alignment.{stk,aln}" in str(exc.value)


def test_workflow_settings_default_and_override_rnaz_no_shuffle() -> None:
    assert WorkflowSettings.from_config({}).rnaz_no_shuffle is True
    assert WorkflowSettings.from_config({"rnaz_no_shuffle": False}).rnaz_no_shuffle is False


def test_verify_run_consistency_accepts_identical_deterministic_outputs(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    for run_dir in (left, right):
        (run_dir / "Lalifold" / "len_150" / "split").mkdir(parents=True)
        (run_dir / "generated_files" / "stk" / "len_150").mkdir(parents=True)
        (run_dir / "Lalifold" / "len_150" / "RC_150_0001.stk").write_text(
            "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n",
            encoding="utf-8",
        )
        (run_dir / "Lalifold" / "len_150" / "split" / "manifest.txt").write_text(
            "RC_150_0001_aln_1_10.stk\n",
            encoding="utf-8",
        )
        (run_dir / "generated_files" / "stk" / "len_150" / "manifest.txt").write_text(
            "RC_150_0001_aln_1_10.stk\n",
            encoding="utf-8",
        )
        split_text = "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n"
        (run_dir / "Lalifold" / "len_150" / "split" / "RC_150_0001_aln_1_10.stk").write_text(
            split_text,
            encoding="utf-8",
        )
        (run_dir / "generated_files" / "stk" / "len_150" / "RC_150_0001_aln_1_10.stk").write_text(
            split_text,
            encoding="utf-8",
        )

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.verify_run_consistency",
            str(left),
            str(right),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0
    assert "Deterministic candidate-generation outputs are identical." in result.stdout


def test_verify_run_consistency_compares_the_real_arm_of_a_calibrated_run(tmp_path: Path) -> None:
    """With the null arm on, every output moves under arms/real/. Comparing the
    run directory itself found no window lengths there, compared nothing, and
    reported the two runs as identical."""
    left, right = tmp_path / "left", tmp_path / "right"
    _write_run_tree(left / "arms" / "real", "ACGU")
    _write_run_tree(right / "arms" / "real", "AGGU")

    result = _run_verify(left, right)

    assert result.returncode == 1
    assert "len_150: cleaned alignments: content differs" in result.stdout


def test_verify_run_consistency_accepts_identical_calibrated_runs(tmp_path: Path) -> None:
    left, right = tmp_path / "left", tmp_path / "right"
    _write_run_tree(left / "arms" / "real", "ACGU")
    _write_run_tree(right / "arms" / "real", "ACGU")

    result = _run_verify(left, right)

    assert result.returncode == 0
    assert "Deterministic candidate-generation outputs are identical." in result.stdout


def test_verify_run_consistency_refuses_to_call_empty_runs_identical(tmp_path: Path) -> None:
    """Having compared nothing is not the same as having found no difference."""
    left, right = tmp_path / "left", tmp_path / "right"
    left.mkdir()
    right.mkdir()

    result = _run_verify(left, right)

    assert result.returncode == 1
    assert "no window lengths to compare" in result.stdout
    assert "identical" not in result.stdout


def test_verify_run_consistency_reports_changed_cleaned_alignment(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    for run_dir in (left, right):
        (run_dir / "Lalifold" / "len_150" / "split").mkdir(parents=True)
        (run_dir / "generated_files" / "stk" / "len_150").mkdir(parents=True)
        (run_dir / "Lalifold" / "len_150" / "RC_150_0001.stk").write_text(
            "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n",
            encoding="utf-8",
        )
        (run_dir / "Lalifold" / "len_150" / "split" / "manifest.txt").write_text(
            "RC_150_0001_aln_1_10.stk\n",
            encoding="utf-8",
        )
        (run_dir / "generated_files" / "stk" / "len_150" / "manifest.txt").write_text(
            "RC_150_0001_aln_1_10.stk\n",
            encoding="utf-8",
        )
        (run_dir / "Lalifold" / "len_150" / "split" / "RC_150_0001_aln_1_10.stk").write_text(
            "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n",
            encoding="utf-8",
        )

    (left / "generated_files" / "stk" / "len_150" / "RC_150_0001_aln_1_10.stk").write_text(
        "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA ACGU\n//\n",
        encoding="utf-8",
    )
    (right / "generated_files" / "stk" / "len_150" / "RC_150_0001_aln_1_10.stk").write_text(
        "# STOCKHOLM 1.0\n#=GF ID RC_150_0001_aln_1_10\nseqA AGGU\n//\n",
        encoding="utf-8",
    )

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.verify_run_consistency",
            str(left),
            str(right),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 1
    assert "len_150: cleaned alignments: content differs" in result.stdout


def test_cli_workflow_smoke_test_with_fake_rnalalifold(tmp_path: Path) -> None:
    input_alignment = tmp_path / "my_input.stk"
    export_dir = tmp_path / "bundle"
    input_alignment.write_text(
        "# STOCKHOLM 1.0\n#=GF ID dummy_input\nseq1 ACGU\nseq2 ACGU\nseq3 ACGU\n//\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_rnalalifold(bin_dir)
    write_fake_postprocess_tools(bin_dir)

    (tmp_path / "config.yaml").write_text(
        textwrap.dedent(
            """\
            lalifold_base: Lalifold
            maxbpspan:
              - 100
              - 200
            lalifold_threads: 1
            do_cm: false
            do_locarnate: false
            # do_png is on here -- and only here -- because true is the shipped
            # default, so the ps2eps -> epstopdf -> magick chain is what a
            # default run actually executes.
            do_png: true
            """
        ),
        encoding="utf-8",
    )

    env = subprocess_env()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["HOME"] = str(tmp_path)
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    env["TMPDIR"] = str(tmp_path / ".tmp")
    (tmp_path / ".cache").mkdir()
    (tmp_path / ".tmp").mkdir()

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.cli",
            "--input-alignment",
            str(input_alignment),
            "--output-dir",
            str(tmp_path),
            "--export-bundle",
            str(export_dir),
            "--maxbpspan",
            "150",
            "--rscape",
            "--cores",
            "1",
            "--",
            "--configfile",
            str(tmp_path / "config.yaml"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "Lalifold" / "len_150" / "RC_150_0001.stk").is_file()
    assert not (tmp_path / "Lalifold" / "len_100").exists()
    assert not (tmp_path / "Lalifold" / "len_200").exists()
    assert (tmp_path / "generated_files" / "remgap" / "len_150" / "manifest.txt").is_file()
    assert (tmp_path / "generated_files" / "strip" / "len_150" / "manifest.txt").is_file()
    assert (tmp_path / "generated_files" / "rscape" / "len_150" / "manifest.txt").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.log").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.log.csv").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.md").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.nr.csv").is_file()
    assert not (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.html").exists()
    refold_out = read_text(
        tmp_path / "generated_files" / "refold" / "len_150" / "RC_150_0001_aln_1_12_refold.out"
    )
    # One record per sequence, as RNAfold -C wrote them before the module took
    # the leg over; the rule no longer needs its own empty-output guard because
    # refold fails on an alignment it cannot read.
    assert refold_out.startswith("> ")
    assert refold_out.count("> ") == 2
    assert _REFOLD_RECORD.search(refold_out), refold_out
    png_dir = tmp_path / "generated_files" / "png" / "len_150"
    assert (png_dir / "manifest.txt").is_file()
    assert (png_dir / "RC_150_0001_aln_1_12_aln.png").is_file()
    assert (png_dir / "RC_150_0001_aln_1_12_ss.png").is_file()
    assert "fake RNALalifold completed for window 150 format S" in read_text(
        tmp_path / "Lalifold" / "len_150" / "RNALalifold.out"
    )
    assert not (tmp_path / "alirna.ps").exists()
    assert not any((tmp_path / "generated_files" / "rnaalifold").glob("**/alirna.ps"))
    assert (
        tmp_path
        / "generated_files"
        / "rnaalifold"
        / "len_150"
        / "RC_150_0001_aln_1_12"
        / "RC_150_0001_aln_1_12_aln.pdf"
    ).is_file()


def test_cli_workflow_smoke_test_accepts_clustal_input(tmp_path: Path) -> None:
    input_alignment = tmp_path / "my_input.aln"
    input_alignment.write_text(
        "CLUSTAL W\n\nseq1 ACGU\nseq2 ACGU\nseq3 ACGU\n***\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_rnalalifold(bin_dir)
    write_fake_postprocess_tools(bin_dir)

    (tmp_path / "config.yaml").write_text(
        textwrap.dedent(
            """\
            lalifold_base: Lalifold
            maxbpspan:
              - 100
              - 200
            lalifold_threads: 1
            do_cm: false
            do_locarnate: false
            do_png: false
            """
        ),
        encoding="utf-8",
    )

    env = subprocess_env()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["HOME"] = str(tmp_path)
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    env["TMPDIR"] = str(tmp_path / ".tmp")
    (tmp_path / ".cache").mkdir()
    (tmp_path / ".tmp").mkdir()

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.cli",
            "--input-alignment",
            str(input_alignment),
            "--output-dir",
            str(tmp_path),
            "--maxbpspan",
            "150",
            "--rscape",
            "--cores",
            "1",
            "--",
            "--configfile",
            str(tmp_path / "config.yaml"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert "fake RNALalifold completed for window 150 format C" in read_text(
        tmp_path / "Lalifold" / "len_150" / "RNALalifold.out"
    )
    assert (tmp_path / "generated_files" / "rscape" / "len_150" / "RC_150_0001_aln_1_12.power").is_file()
    assert "RC_150_0001_aln_1_12" in read_text(
        tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.log"
    )
    assert "rscape_covary_count" in read_text(
        tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.log.csv"
    )


def test_candidate_paths_prefixed_by_arm_directory() -> None:
    plain = CandidatePaths(wlen=100, file="RC_100_0001_aln_1_12")
    armed = CandidatePaths(wlen=100, file="RC_100_0001_aln_1_12", arm_prefix="arms/null_000/")

    assert plain.rnaz_json == "generated_files/rnaz/len_100/RC_100_0001_aln_1_12.rnaz.json"
    assert armed.rnaz_json == "arms/null_000/" + plain.rnaz_json
    assert armed.aln_pdf == "arms/null_000/" + plain.aln_pdf
    assert armed.split == "arms/null_000/" + plain.split


def test_snakefile_defines_every_rule_exactly_once() -> None:
    import re as _re

    text = read_text(Path("snakefile"))
    names = _re.findall(r"^\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", text, _re.MULTILINE)
    duplicates = sorted({name for name in names if names.count(name) > 1})
    assert not duplicates, f"duplicated rule definitions: {duplicates}"

    # No parallel arm-specific rules: the arms must traverse the same DAG.
    assert not [name for name in names if name.endswith(("_null", "_real"))]
    assert "make_arm_alignment" in names
    assert "calibrate" in names


def test_snakefile_has_no_duplicated_rule_bodies() -> None:
    import re as _re

    text = read_text(Path("snakefile"))
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.MULTILINE)
    bodies: dict[str, str] = {}
    for block in blocks:
        match = _re.match(r"\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", block)
        if not match:
            continue
        body = "\n".join(line.strip() for line in block.split("\n", 1)[1].splitlines() if line.strip())
        if not body:
            continue
        for other, other_body in bodies.items():
            assert body != other_body, f"rules {match.group(1)} and {other} have identical bodies"
        bodies[match.group(1)] = body


def test_config_null_section_is_quoted_and_disabled_by_default() -> None:
    import yaml as _yaml

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    # An unquoted null: key would land under None and silently do nothing.
    assert None not in payload
    assert payload["null"]["method"] == "none"
    assert payload["calibration"]["stage1_rnaz_prob"] <= payload["calibration"]["rnaz_prob_threshold"]


def test_cli_workflow_null_arm_smoke_test(tmp_path: Path) -> None:
    """Whole DAG, including two null replicates, on a toy set.

    This is the CI target from the calibration handoff: it proves the real and
    null arms traverse the same rules and that funnel.tsv gets one row per
    filter stage per arm.
    """
    if shutil.which("perl") is None:  # pragma: no cover
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    input_alignment = tmp_path / "my_input.stk"
    input_alignment.write_text(TOY_ALIGNMENT, encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_rnalalifold(bin_dir)
    write_fake_postprocess_tools(bin_dir)
    randomize = write_fake_randomize_aln(bin_dir)

    # The fake RNALalifold reports one window over columns 1-12 in every arm,
    # so a truth element there is recovered by the real arm -- and by the null
    # arms too, which is exactly what the baseline exists to expose.
    truth = tmp_path / "truth.tsv"
    truth.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "toy_element\txrRNA\ttoy\t1\t12\tsynthetic\n",
        encoding="utf-8",
    )

    (tmp_path / "config.yaml").write_text(
        textwrap.dedent(
            f"""\
            lalifold_base: Lalifold
            maxbpspan:
              - 150
            lalifold_threads: 1
            do_cm: false
            do_locarnate: false
            do_png: false
            "null":
              method: rnazRandomizeAln
              replicates: 2
              seed: 20261101
              two_stage: true
            calibration:
              rnaz_prob_threshold: 0.9
              alifoldz_threshold: -2.0
              rscape_min_pairs: 1
              stage1_rnaz_prob: 0.5
            benchmark_truth: {truth}
            tools:
              rnaz_randomize_aln: {randomize}
            """
        ),
        encoding="utf-8",
    )

    env = subprocess_env()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["HOME"] = str(tmp_path)
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    env["TMPDIR"] = str(tmp_path / ".tmp")
    (tmp_path / ".cache").mkdir()
    (tmp_path / ".tmp").mkdir()

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.cli",
            "--input-alignment",
            str(input_alignment),
            "--output-dir",
            str(tmp_path),
            "--maxbpspan",
            "150",
            "--benchmark",
            "--cores",
            "1",
            "--",
            "--configfile",
            str(tmp_path / "config.yaml"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    # Every arm traversed the same rules, under its own prefix.
    for arm in ["real", "null_000", "null_001"]:
        assert (tmp_path / "arms" / arm / "input_alignment.stk").is_file()
        assert (tmp_path / "arms" / arm / "Lalifold" / "len_150" / "RC_150_0001.stk").is_file()
        assert (
            tmp_path / "arms" / arm / "generated_files" / "summary" / "len_150" / "RNAConSnake.log.csv"
        ).is_file()
    assert not (tmp_path / "Lalifold").exists()

    # The null replicates differ from the real alignment and from each other.
    real_arm = read_text(tmp_path / "arms" / "real" / "input_alignment.stk")
    null_0 = read_text(tmp_path / "arms" / "null_000" / "input_alignment.stk")
    null_1 = read_text(tmp_path / "arms" / "null_001" / "input_alignment.stk")
    assert real_arm == TOY_ALIGNMENT
    assert null_0 != real_arm
    assert null_1 != real_arm
    assert null_0 != null_1

    calibration = tmp_path / "results" / "calibration"
    funnel = read_text(calibration / "funnel.tsv")
    for arm in ["real", "null_000", "null_001"]:
        for stage in ["windows", "loci", "rnaz", "alifoldz", "rscape", "cascade"]:
            assert f"\t{arm}\t150\t{stage}\t" in funnel
    assert "# fdr_conditional_on_stage_one\ttrue" in funnel

    assert (calibration / "qvalues.tsv").is_file()
    assert (calibration / "score_dists.tsv").is_file()
    summary = json.loads(read_text(calibration / "summary.json"))
    assert summary["null"]["method"] == "rnazRandomizeAln"
    assert summary["null"]["replicates"] == 2
    assert summary["fdr_conditional_on_stage_one"] is True
    assert summary["counting_unit"] == "merged_loci"

    # Reproducibility bookkeeping travels with the calibrated run. This also
    # guards --benchmark being *additive*: naming a target replaces the default
    # one, so asking for the recovery table must not drop the run's outputs.
    assert (tmp_path / "results" / "versions.yaml").is_file()
    assert (tmp_path / "null_pool" / "pool.json").is_file()

    # The positive control ran, and it carries its null baseline: a recovery
    # count without one says nothing about whether the screen found structure.
    recovery = read_text(tmp_path / "results" / "benchmark" / "flavivirus_recovery.tsv")
    assert "# recovered\t1" in recovery
    assert "over 2 null arms" in recovery
    assert "not evidence of detection" in recovery  # both arms recover it here
    assert "toy_element" in recovery


def test_ci_test_profile_builds_the_calibrated_dag(tmp_path: Path) -> None:
    """profiles/test documents itself as the CI / toy run, so it has to keep
    producing the DAG it claims: both arms plus calibration, PNGs off. A dry
    run needs no external tools, so this holds wherever the suite runs."""
    import re as _re

    alignment = tmp_path / "toy.stk"
    alignment.write_text(TOY_ALIGNMENT, encoding="utf-8")
    # input_alignment travels in a configfile: a command-line --config would
    # replace the profile's own config entries wholesale.
    configfile = tmp_path / "config.yaml"
    configfile.write_text(f"input_alignment: {alignment}\n", encoding="utf-8")

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "snakemake",
            "-n",
            "--snakefile",
            str(ROOT / "snakefile"),
            "--profile",
            str(ROOT / "profiles" / "test"),
            "--directory",
            str(tmp_path),
            "--configfile",
            str(configfile),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stdout + result.stderr
    output = result.stdout + result.stderr
    # The real arm and both null replicates, and the step that needs them.
    assert _re.search(r"^make_arm_alignment\s+3$", output, _re.MULTILINE), output
    assert _re.search(r"^calibrate\s+1$", output, _re.MULTILINE), output
    # do_png=False in the profile, so the image branch stays out of the DAG.
    assert "render_pngs_file" not in output


def test_cli_workflow_none_method_keeps_legacy_layout(tmp_path: Path) -> None:
    """`method: none` must be as inert as not having the null arm at all."""
    input_alignment = tmp_path / "my_input.stk"
    input_alignment.write_text(TOY_ALIGNMENT, encoding="utf-8")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_rnalalifold(bin_dir)
    write_fake_postprocess_tools(bin_dir)

    (tmp_path / "config.yaml").write_text(
        textwrap.dedent(
            """\
            lalifold_base: Lalifold
            maxbpspan:
              - 150
            lalifold_threads: 1
            do_cm: false
            do_locarnate: false
            do_png: false
            "null":
              method: none
              replicates: 10
            """
        ),
        encoding="utf-8",
    )

    env = subprocess_env()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["HOME"] = str(tmp_path)
    env["XDG_CACHE_HOME"] = str(tmp_path / ".cache")
    env["TMPDIR"] = str(tmp_path / ".tmp")
    (tmp_path / ".cache").mkdir()
    (tmp_path / ".tmp").mkdir()

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.cli",
            "--input-alignment",
            str(input_alignment),
            "--output-dir",
            str(tmp_path),
            "--maxbpspan",
            "150",
            "--cores",
            "1",
            "--",
            "--configfile",
            str(tmp_path / "config.yaml"),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stdout + result.stderr

    assert (tmp_path / "Lalifold" / "len_150" / "RC_150_0001.stk").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.log.csv").is_file()
    assert not (tmp_path / "arms").exists()
    assert not (tmp_path / "results").exists()
    assert not (tmp_path / "null_pool").exists()


def test_run_checked_never_inherits_stdin(tmp_path: Path) -> None:
    """A tool that reads stdin must hit EOF, not block on the user's terminal.

    Perl's <> falls through to STDIN once @ARGV is exhausted -- alifoldz.pl
    still does. Inheriting an interactive terminal there parks the job forever
    and holds a scheduler slot.
    """
    from rnaconsnake.workflow_helpers import run_checked

    reader = tmp_path / "reader.py"
    reader.write_text(
        "import sys\nsys.stdout.write(repr(sys.stdin.read()))\n",
        encoding="utf-8",
    )
    output = tmp_path / "out.txt"
    run_checked([PYTHON, str(reader)], stdout_path=output)
    assert read_text(output) == "''"


def test_run_checked_still_honours_an_explicit_stdin_path(tmp_path: Path) -> None:
    from rnaconsnake.workflow_helpers import run_checked

    reader = tmp_path / "reader.py"
    reader.write_text("import sys\nsys.stdout.write(sys.stdin.read())\n", encoding="utf-8")
    source = tmp_path / "in.txt"
    source.write_text("payload\n", encoding="utf-8")
    output = tmp_path / "out.txt"
    run_checked([PYTHON, str(reader)], stdin_path=source, stdout_path=output)
    assert read_text(output) == "payload\n"


def test_workflow_subprocesses_do_not_inherit_stdin() -> None:
    """Every subprocess launched by the workflow pins stdin explicitly."""
    import re as _re

    text = read_text(Path("snakefile"))
    for match in _re.finditer(r"subprocess\.(run|Popen)\(", text):
        call = text[match.start() : match.start() + 600]
        assert "stdin=" in call.split("        )")[0], (
            f"subprocess call at offset {match.start()} does not pin stdin"
        )


def test_configurable_rules_declare_their_settings_as_params() -> None:
    """Config read from globals inside `run:` is invisible to Snakemake.

    Without a `params:` declaration the params rerun-trigger cannot see a
    changed setting, so re-running with a different method silently reuses the
    previous output.
    """
    import re as _re

    text = read_text(Path("snakefile"))
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.MULTILINE)
    by_name = {}
    for block in blocks:
        match = _re.match(r"\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", block)
        if match:
            by_name[match.group(1)] = block

    for rule in ["dereplicate_summary", "calibrate", "benchmark_recovery"]:
        assert rule in by_name, f"{rule} is missing from the snakefile"
        assert _re.search(r"^\s+params:", by_name[rule], _re.MULTILINE), (
            f"rule {rule} reads configurable settings but declares no params:, "
            "so changing them will not trigger a re-run"
        )

    # ...and the run body must consume params, not the globals directly,
    # otherwise the declaration is decorative.
    body = by_name["dereplicate_summary"].split("\n    run:", 1)[1]
    assert "params.method" in body
    assert "DEREPLICATE_METHOD" not in body
    assert "PAIR_CONTAINMENT" not in body


def test_versions_report_records_python_and_external_tools(tmp_path: Path) -> None:
    from rnaconsnake.tools.versions import write_versions

    output = tmp_path / "versions.yaml"
    payload = write_versions(output, {"rnalalifold": "definitely-not-installed"})
    assert payload["external_tools"]["rnalalifold"]["path"] == "MISSING"
    text = read_text(output)
    assert text.startswith("# RNAcs toolchain versions")
    assert "python:" in text
    assert "external_tools:" in text
    # Refolding runs through the ViennaRNA Python module, so its version is
    # part of the toolchain a calibrated run has to be reproducible against.
    assert "viennarna_bindings:" in text
    assert payload["viennarna_bindings"]["module"] == "RNA"
    assert payload["viennarna_bindings"]["matches_binaries"] in {"yes", "no", "unknown"}


def test_versions_report_is_valid_yaml_even_with_awkward_tool_paths(tmp_path: Path) -> None:
    """It is written with yaml.safe_dump now. The hand-rolled writer quoted every
    value with `"` and escaped nothing but `"`, so a backslash in a path emitted
    a file that would not parse back."""
    import yaml as _yaml

    from rnaconsnake.tools.versions import write_versions

    output = tmp_path / "versions.yaml"
    awkward = r"C:\tools\RNAlalifold --flag 'x' \"quoted\""
    write_versions(output, {"rnalalifold": awkward})

    parsed = _yaml.safe_load(read_text(output))
    assert parsed["external_tools"]["rnalalifold"]["command"] == awkward
    assert parsed["python"]
    assert parsed["viennarna_bindings"]["module"] == "RNA"


def test_versions_report_states_whether_the_toolchain_is_one_build(tmp_path: Path) -> None:
    """The finished run's provenance answers the question on its own."""
    from rnaconsnake.tools.versions import viennarna_consistency

    bindings = {"version": "2.7.2"}
    matching = {
        "rnalalifold": {"version": "RNALalifold 2.7.2"},
        "rnaalifold": {"version": "RNAalifold 2.7.2"},
    }
    mixed = {"rnalalifold": {"version": "RNALalifold 2.7.2"}, "rnaalifold": {"version": "RNAalifold 2.6.4"}}
    assert viennarna_consistency(bindings, matching) == "yes"
    assert viennarna_consistency(bindings, mixed) == "no"
    assert viennarna_consistency({"version": "unknown"}, matching) == "unknown"
    assert viennarna_consistency(bindings, {"rnalalifold": {"version": "unknown"}}) == "unknown"


def test_version_probes_do_not_write_into_the_working_directory(tmp_path: Path) -> None:
    """ViennaRNA tools drop side outputs wherever they are run, and the version
    probe runs in the user's own directory."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    litterbug = bin_dir / "litterbug"
    litterbug.write_text(
        "#!/usr/bin/env python3\n"
        "from pathlib import Path\n"
        "Path('side_output.ps').write_text('%!PS\\n', encoding='utf-8')\n"
        "print('litterbug 1.0')\n",
        encoding="utf-8",
    )
    litterbug.chmod(0o755)

    here = tmp_path / "cwd"
    here.mkdir()
    monkeypatched_path = f"{bin_dir}:{os.environ['PATH']}"
    original_path, original_cwd = os.environ["PATH"], os.getcwd()
    os.environ["PATH"] = monkeypatched_path
    os.chdir(here)
    try:
        assert cli.probe_version("litterbug") == "litterbug 1.0"
        from rnaconsnake.tools.versions import probe

        assert probe("litterbug", ["--version"])["version"] == "litterbug 1.0"
    finally:
        os.chdir(original_cwd)
        os.environ["PATH"] = original_path

    assert not list(here.iterdir()), "a version probe wrote into the working directory"


def test_the_refold_leg_no_longer_needs_perl_or_rnafold() -> None:
    """The refold leg is ours now. refold.pl and the RNAfold binary must not
    creep back in through the workflow, the config, or the preflight."""
    import yaml as _yaml

    # The rule may still name the pipe it replaced; it must not invoke it.
    snakefile = read_text(Path("snakefile"))
    assert 'command_tokens("refold", "python3 -m rnaconsnake.tools.refold")' in snakefile
    assert 'command_tokens("rnafold"' not in snakefile
    assert '"refold.pl"' not in snakefile

    config = _yaml.safe_load(read_text(Path("config.yaml")))
    assert config["tools"]["refold"] == "python3 -m rnaconsnake.tools.refold"
    assert "rnafold" not in config["tools"]

    # Not an executable to look for on PATH any more.
    assert "refold" not in cli.DEFAULT_RUNTIME_TOOLS
    assert "rnafold" not in cli.DEFAULT_RUNTIME_TOOLS
    assert cli.DEFAULT_TOOL_COMMANDS["refold"].startswith("python3 -m rnaconsnake")
