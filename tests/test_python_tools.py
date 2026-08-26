from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from rnaconsnake import cli
from rnaconsnake.workflow_helpers import CandidatePaths, WorkflowSettings, initial_alignment_format_code

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
PYTHON = sys.executable


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    pythonpath = str(ROOT / "src")
    env["PYTHONPATH"] = pythonpath + (f":{env['PYTHONPATH']}" if env.get("PYTHONPATH") else "")
    return env


def write_fake_rnalalifold(bin_dir: Path) -> None:
    script = bin_dir / "RNALalifold"
    script.write_text(
        textwrap.dedent(
            """\
            #!/usr/bin/env python3
            import sys
            from pathlib import Path

            args = sys.argv[1:]
            if "--version" in args:
                import RNA

                print("RNALalifold " + RNA.__version__)
                raise SystemExit(0)
            prefix = "RC"
            wlen = "100"
            input_format = ""
            for i, arg in enumerate(args):
                if arg == "--id-prefix" and i + 1 < len(args):
                    prefix = args[i + 1]
                if arg == "-L" and i + 1 < len(args):
                    wlen = args[i + 1]
                if arg == "-f" and i + 1 < len(args):
                    input_format = args[i + 1]

            out_path = Path.cwd() / f"{prefix}_0001.stk"
            out_path.write_text(
                "\\n".join(
                    [
                        "# STOCKHOLM 1.0",
                        f"#=GF ID {prefix}_0001_aln_1_12",
                        "seqA ACGUACGU----",
                        "seqB ACGUACGU----",
                        "#=GC SS_cons <<<<....>>>>",
                        "//",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            sys.stdout.write(f"fake RNALalifold completed for window {wlen} format {input_format}\\n")
            """
        ),
        encoding="utf-8",
    )
    script.chmod(0o755)


def write_fake_postprocess_tools(bin_dir: Path) -> None:
    tools = {
        "esl-reformat": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[-1])
print("CLUSTAL W")
print()
for line in inp.read_text(encoding="utf-8").splitlines():
    if not line or line.startswith("#") or line == "//":
        continue
    name, seq = line.split(None, 1)
    print(f"{name} {seq}")
print("***")
""",
        "RNAz": """\
#!/usr/bin/env python3
print("Mean z-score: -3.21")
print("Mean MFE: -14.10")
print("Structure conservation index: 0.58")
print("SVM RNA-class probability: 0.95")
""",
        "alifoldz.pl": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

# The real alifoldz.pl shells out to RNAalifold without --noPS, so an
# undeclared "alirna.ps" lands in whatever directory it was run from. Emit one
# too, so the workflow's confinement of that stray file is actually exercised.
(Path.cwd() / "alirna.ps").write_text("%!PS stray from alifoldz\\n", encoding="utf-8")
_ = sys.stdin.read()
print("#           Input: 3 sequences of 41 columns")
print("fake alifoldz header")
print("-3.21")
""",
        "RNAalifold": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if "--version" in args:
    # The real binary prints and exits; it does not fold, and must not write.
    # The version is the module's: RNAcs requires one ViennaRNA build.
    import RNA
    print("RNAalifold " + RNA.__version__)
    raise SystemExit(0)
prefix = "fake"
for i, arg in enumerate(args):
    if arg == "--id-prefix" and i + 1 < len(args):
        prefix = args[i + 1]

cwd = Path.cwd()
(cwd / "RNAalifold_results.stk").write_text(
    "\\n".join(
        [
            "# STOCKHOLM 1.0",
            f"#=GF ID {prefix}",
            "seqA ACGUACGU----",
            "seqB ACGUACGU----",
            "#=GC SS_cons <<<<....>>>>",
            "//",
            "",
        ]
    ),
    encoding="utf-8",
)
(cwd / f"{prefix}_0001_ali.out").write_text("2 a b c d 12\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_dp.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_aln.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / f"{prefix}_0001_ss.ps").write_text("%!PS\\n", encoding="utf-8")
(cwd / "alirna.ps").write_text("%!PS stray\\n", encoding="utf-8")
print("fake RNAalifold run")
""",
        "ps2eps": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = inp.with_suffix(".eps")
out.write_text("%!EPS\\n", encoding="utf-8")
""",
        "epstopdf": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = inp.with_suffix(".pdf")
out.write_text("%PDF-FAKE\\n", encoding="utf-8")
""",
        "magick": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

out = Path(sys.argv[-1])
out.write_text("PNG-FAKE\\n", encoding="utf-8")
""",
        "R-scape": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

inp = Path(sys.argv[-1])
(Path.cwd() / f"{inp.stem}.power").write_text("# BPAIRS observed to covary 1\\n", encoding="utf-8")
""",
    }

    for name, content in tools.items():
        script = bin_dir / name
        script.write_text(textwrap.dedent(content), encoding="utf-8")
        script.chmod(0o755)


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
    assert paths.refold_json == "generated_files/refold/len_150/RC_150_0001_aln_1_12.refold.json"
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
    with pytest.raises(ValueError) as exc:
        initial_alignment_format_code(None)
    assert "input_alignment.{stk,aln}" in str(exc.value)


def test_workflow_settings_default_and_override_rnaz_no_shuffle() -> None:
    assert WorkflowSettings.from_config({}).rnaz_no_shuffle is True
    assert WorkflowSettings.from_config({"rnaz_no_shuffle": False}).rnaz_no_shuffle is False


def test_export_bundle_writes_spec_files_and_copies_artifacts(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    candidate_id = "RC_150_0001_aln_1_12"
    paths = CandidatePaths(wlen=150, file=candidate_id)

    summary_dir = run_dir / "generated_files" / "summary" / "len_150"
    summary_dir.mkdir(parents=True)
    (summary_dir / f"{candidate_id}.summary.json").write_text(
        json.dumps(
            {
                "wbn": candidate_id,
                "nrseq": "2",
                "alilen": "12",
                "maxcovarval": "5",
                "maxcovarcount": "1",
                "rscape_covary_count": "1",
                "rnazprob": "0.95",
                "sci": "0.58",
                "consensus_mfe": "-14.10",
                "alifoldzscore": "-3.21",
                "alifold_consstruc": "<<>>",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (summary_dir / "RNAConSnake.log").write_text("log\n", encoding="utf-8")
    (summary_dir / "RNAConSnake.log.csv").write_text("header\n", encoding="utf-8")
    (summary_dir / "RNAConSnake.md").write_text("# summary\n", encoding="utf-8")

    for relpath, content in {
        paths.split: "# STOCKHOLM 1.0\n//\n",
        paths.stk: "# STOCKHOLM 1.0\n//\n",
        paths.aln: "CLUSTAL W\n",
        paths.refold_json: '{"alifold_consstruc": "<<>>"}\n',
        paths.aln_pdf: "%PDF-FAKE\n",
        paths.ss_pdf: "%PDF-FAKE\n",
        paths.rscape_power: "# BPAIRS observed to covary 1\n",
        paths.rscape_sto_pdf: "%PDF-FAKE\n",
    }.items():
        target = run_dir / relpath
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8")

    cm_status = run_dir / paths.cm_status_json
    cm_status.parent.mkdir(parents=True, exist_ok=True)
    cm_status.write_text(
        '{"built": true, "cm": "generated_files/cm/len_150/RC_150_0001_aln_1_12.cm"}\n', encoding="utf-8"
    )
    (run_dir / "generated_files" / "cm" / "len_150" / f"{candidate_id}.cm").write_text(
        "CM\n", encoding="utf-8"
    )

    bundle_dir = tmp_path / "bundle"
    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.export_bundle",
            str(run_dir),
            str(bundle_dir),
            "--input-alignment",
            str(tmp_path / "np.mlocarna2020.stk"),
            "--feature-type",
            "other",
            "--overwrite",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stderr
    manifest = json.loads(read_text(bundle_dir / "manifest.json"))
    assert manifest["project"] == "RNAConSnake"
    # 1.1.0 added the de-replication and calibration columns and manifest blocks.
    assert manifest["export_schema_version"] == "1.1.0"
    assert manifest["dereplication"]["applied"] is False
    assert manifest["calibration"]["applied"] is False
    assert manifest["features_file"] == "features.csv"
    assert manifest["candidates_file"] == "candidates.csv"
    assert manifest["artifacts_file"] == "artifacts.csv"

    features_text = read_text(bundle_dir / "features.csv")
    candidates_text = read_text(bundle_dir / "candidates.csv")
    artifacts_text = read_text(bundle_dir / "artifacts.csv")

    assert "feature_id,feature_label,feature_type,candidate_count,passing_count" in features_text
    assert "candidate_id,feature_id,start,end,length,n_sequences,status,score_primary" in candidates_text
    assert "rnazprob" in candidates_text
    assert "sci" in candidates_text
    assert "consensus_mfe" in candidates_text
    assert "alifoldzscore" in candidates_text
    assert (
        "artifact_scope,artifact_owner_id,artifact_type,artifact_label,file_format,path,is_optional"
        in artifacts_text
    )
    assert "np.mlocarna2020" in features_text
    assert candidate_id in candidates_text
    assert "0.95" in candidates_text
    assert "0.58" in candidates_text
    assert "-14.10" in candidates_text
    assert "-3.21" in candidates_text
    assert "alignment_plot_pdf" in artifacts_text
    assert "consensus_plot_pdf" in artifacts_text
    assert "rscape_pdf" in artifacts_text
    assert "covariance_model" in artifacts_text
    assert (bundle_dir / "methods.md").is_file()
    assert (bundle_dir / "feature_summaries" / "np.mlocarna2020.md").is_file()
    assert (bundle_dir / "candidate_summaries" / f"{candidate_id}.md").is_file()
    assert (
        bundle_dir / "files" / "np.mlocarna2020" / candidate_id / f"{candidate_id}.alignment_plot.pdf"
    ).is_file()
    assert (
        bundle_dir / "files" / "np.mlocarna2020" / candidate_id / f"{candidate_id}.consensus_plot.pdf"
    ).is_file()


def test_render_reports_writes_log_csv_and_markdown(tmp_path: Path) -> None:
    summary = tmp_path / "one.summary.json"
    summary.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec1",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "1",
              "maxcovarcount": "2",
              "rscape_covary_count": "1",
              "rnazprob": "0.95",
              "alifoldzscore": "-3.21",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "write-summary-outputs",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "RNAConSnake.log"),
            "--csv",
            str(tmp_path / "RNAConSnake.log.csv"),
            "--markdown",
            str(tmp_path / "RNAConSnake.md"),
            str(summary),
        ],
        check=True,
        env=subprocess_env(),
    )

    assert "rec1" in read_text(tmp_path / "RNAConSnake.log")
    assert "wbn,nrseq,alilen" in read_text(tmp_path / "RNAConSnake.log.csv")
    assert "# RNAConSnake Summary: len_100" in read_text(tmp_path / "RNAConSnake.md")


def test_extract_rnaz_reads_probability_sci_and_consensus_mfe(tmp_path: Path) -> None:
    rnaz = tmp_path / "one.rnaz.txt"
    out = tmp_path / "one.rnaz.json"
    rnaz.write_text(
        "\n".join(
            [
                "Mean MFE: -14.10",
                "Structure conservation index: 0.58",
                "SVM RNA-class probability: 0.95",
                "",
            ]
        ),
        encoding="utf-8",
    )

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "extract-rnaz",
            "--input",
            str(rnaz),
            "--output",
            str(out),
        ],
        check=True,
        env=subprocess_env(),
    )

    payload = json.loads(read_text(out))
    assert payload == {"consensus_mfe": "-14.10", "rnazprob": "0.95", "sci": "0.58"}


def test_write_summary_outputs_sorts_records_by_maxcovar_then_alifoldz(tmp_path: Path) -> None:
    """Ties on covariation break towards the *more negative* AlifoldZ, which is
    the more significant one -- the ordering de-replication also uses."""
    summary_a = tmp_path / "a.summary.json"
    summary_b = tmp_path / "b.summary.json"
    summary_c = tmp_path / "c.summary.json"

    summary_a.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec_a",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "3",
              "maxcovarcount": "1",
              "rscape_covary_count": "",
              "rnazprob": "0.5",
              "alifoldzscore": "-1.0",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )
    summary_b.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec_b",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "5",
              "maxcovarcount": "1",
              "rscape_covary_count": "",
              "rnazprob": "0.5",
              "alifoldzscore": "-2.0",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )
    summary_c.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec_c",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "5",
              "maxcovarcount": "1",
              "rscape_covary_count": "",
              "rnazprob": "0.5",
              "alifoldzscore": "-3.0",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "write-summary-outputs",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "RNAConSnake.log"),
            "--csv",
            str(tmp_path / "RNAConSnake.log.csv"),
            "--markdown",
            str(tmp_path / "RNAConSnake.md"),
            str(summary_a),
            str(summary_b),
            str(summary_c),
        ],
        check=True,
        env=subprocess_env(),
    )

    csv_lines = read_text(tmp_path / "RNAConSnake.log.csv").splitlines()
    # rec_b and rec_c tie on maxcovarval 5; rec_c's AlifoldZ of -3.0 beats -2.0.
    assert csv_lines[1].startswith("rec_c,")
    assert csv_lines[2].startswith("rec_b,")
    assert csv_lines[3].startswith("rec_a,")


def test_alifoldz_ranks_the_same_way_everywhere(tmp_path: Path) -> None:
    """A more negative AlifoldZ is the more significant one. Every place that
    orders candidates by it has to agree, or the export bundle and the run's
    own reports disagree with the locus representative de-replication picked."""
    from rnaconsnake.export_bundle import SummaryRecord, sort_summary_records
    from rnaconsnake.tools.dereplicate import rank_key
    from rnaconsnake.tools.legacy_postprocess import _report_sort_key

    weak = {"maxcovarval": "5", "alifoldzscore": "-2.0", "rnazprob": "0.9"}
    strong = {"maxcovarval": "5", "alifoldzscore": "-3.0", "rnazprob": "0.9"}
    missing = {"maxcovarval": "5", "alifoldzscore": "NA", "rnazprob": "0.9"}

    assert _report_sort_key(strong) > _report_sort_key(weak) > _report_sort_key(missing)
    assert rank_key(strong) > rank_key(weak) > rank_key(missing)

    records = [
        SummaryRecord(wlen=100, values=dict(values, wbn=name), summary_path=tmp_path / name)
        for name, values in [("weak", weak), ("strong", strong), ("missing", missing)]
    ]
    assert [record.candidate_id for record in sort_summary_records(records)] == [
        "strong",
        "weak",
        "missing",
    ]


def test_legacy_postprocess_removed_render_reports_subcommand(tmp_path: Path) -> None:
    summary = tmp_path / "one.summary.json"
    summary.write_text('{"wbn":"rec1"}\n', encoding="utf-8")

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "render-reports",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "RNAConSnake.log"),
            "--csv",
            str(tmp_path / "RNAConSnake.log.csv"),
            "--markdown",
            str(tmp_path / "RNAConSnake.md"),
            str(summary),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode != 0
    assert "invalid choice" in result.stderr
    assert "write-summary-outputs" in result.stderr


def test_split_stockholm_matches_expected_outputs(tmp_path: Path) -> None:
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.split_stockholm",
            "-a",
            str(FIXTURES / "multi_input.stk"),
        ],
        cwd=tmp_path,
        check=True,
        env=subprocess_env(),
    )

    assert read_text(tmp_path / "rec1.stk") == read_text(FIXTURES / "expected_split_rec1.stk")
    assert read_text(tmp_path / "acc2.stk") == read_text(FIXTURES / "expected_split_acc2.stk")


def test_split_stockholm_matches_real_lalifold_sample(tmp_path: Path) -> None:
    real = FIXTURES / "real"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.split_stockholm",
            "-a",
            str(real / "lalifold_multi_sample.stk"),
        ],
        cwd=tmp_path,
        check=True,
        env=subprocess_env(),
    )

    assert read_text(tmp_path / "RC_100_0001_aln_2294_2375.stk") == read_text(
        real / "RC_100_0001_aln_2294_2375.stk"
    )
    assert read_text(tmp_path / "RC_100_0001_aln_4042_4067.stk") == read_text(
        real / "RC_100_0001_aln_4042_4067.stk"
    )


def test_remove_gaponly_matches_expected_output(tmp_path: Path) -> None:
    out = tmp_path / "out.stk"
    with out.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.remove_gaponly",
                "-a",
                str(FIXTURES / "expected_split_rec1.stk"),
                "-i",
                "stockholm",
                "-r",
                "0.5",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    assert read_text(out) == read_text(FIXTURES / "expected_remgap_rec1.stk")


def test_remove_gaponly_filters_sequences_by_max_n(tmp_path: Path) -> None:
    inp = tmp_path / "with_n.stk"
    inp.write_text(
        textwrap.dedent(
            """\
            # STOCKHOLM 1.0
            seqA ACGUACGU
            seqB ACNUNCGU
            seqC ACGUACGU
            //
            """
        ),
        encoding="utf-8",
    )

    out_zero = tmp_path / "out_n0.stk"
    with out_zero.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.remove_gaponly",
                "-a",
                str(inp),
                "-i",
                "stockholm",
                "-r",
                "1.0",
                "-n",
                "0",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    out_two = tmp_path / "out_n2.stk"
    with out_two.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.remove_gaponly",
                "-a",
                str(inp),
                "-i",
                "stockholm",
                "-r",
                "1.0",
                "-n",
                "2",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    zero_text = read_text(out_zero)
    two_text = read_text(out_two)
    assert "seqB" not in zero_text
    assert "seqB" in two_text


def test_real_lalifold_sample_matches_expected_remgap_and_strip(tmp_path: Path) -> None:
    real = FIXTURES / "real"
    remgap_out = tmp_path / "remgap.stk"
    strip_out = tmp_path / "strip.stk"

    with remgap_out.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.remove_gaponly",
                "-a",
                str(real / "RC_100_0001_aln_2294_2375.stk"),
                "-i",
                "stockholm",
                "-r",
                "0.5",
                "-n",
                "999",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    with strip_out.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.strip_aln",
                "-a",
                str(remgap_out),
                "-f",
                "S",
                "--nosingle",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    assert read_text(remgap_out) == read_text(real / "expected_remgap_2294_2375.stk")
    assert read_text(strip_out) == read_text(real / "expected_strip_2294_2375.stk")


def test_strip_aln_nosingle_preserves_all_identical_sequences(tmp_path: Path) -> None:
    out = tmp_path / "out.stk"
    with out.open("w", encoding="utf-8") as handle:
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.strip_aln",
                "-a",
                str(FIXTURES / "expected_split_acc2.stk"),
                "-f",
                "S",
                "--nosingle",
            ],
            check=True,
            stdout=handle,
            env=subprocess_env(),
        )

    assert read_text(out) == read_text(FIXTURES / "expected_strip_acc2_nosingle.stk")


def test_alifold_maxcovar_matches_legacy_behavior() -> None:
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.tools.alifold_maxcovar"],
        input="  1 x 2 a b c d e\n  2 x 1 a b c d --:4\n",
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert result.returncode == 3
    assert result.stdout == "3\tgapali_1"


def test_extract_rscape_reads_covary_count(tmp_path: Path) -> None:
    power = tmp_path / "one.power"
    out = tmp_path / "one.rscape.json"
    power.write_text("# BPAIRS observed to covary 1\n", encoding="utf-8")

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "extract-rscape",
            "--input",
            str(power),
            "--output",
            str(out),
        ],
        check=True,
        env=subprocess_env(),
    )

    assert '"rscape_covary_count": "1"' in read_text(out)


def test_extract_rscape_reads_zero_covary_count(tmp_path: Path) -> None:
    power = tmp_path / "zero.power"
    out = tmp_path / "zero.rscape.json"
    power.write_text("# BPAIRS observed to covary 0\n", encoding="utf-8")

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "extract-rscape",
            "--input",
            str(power),
            "--output",
            str(out),
        ],
        check=True,
        env=subprocess_env(),
    )

    assert '"rscape_covary_count": "0"' in read_text(out)


def test_extract_rscape_returns_empty_without_covary_line(tmp_path: Path) -> None:
    power = tmp_path / "one.power"
    out = tmp_path / "one.rscape.json"
    power.write_text("# no covary summary here\n", encoding="utf-8")

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "extract-rscape",
            "--input",
            str(power),
            "--output",
            str(out),
        ],
        check=True,
        env=subprocess_env(),
    )

    assert '"rscape_covary_count": ""' in read_text(out)


def test_failed_alifoldz_fallback_is_non_numeric(tmp_path: Path) -> None:
    """A failed alifoldz must not look like a real, unremarkable z-score.

    The snakefile writes "NA" rather than "0.0" so a crashed run cannot enter
    the summary tables or the FDR calibration as a plausible score.
    """
    snakefile = read_text(Path("snakefile"))
    assert '{"alifoldzscore": "NA"}' in snakefile
    assert '{"alifoldzscore": "0.0"}' not in snakefile


def test_write_summary_outputs_preserves_rscape_zero_and_na(tmp_path: Path) -> None:
    zero_summary = tmp_path / "zero.summary.json"
    na_summary = tmp_path / "na.summary.json"

    zero_summary.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec_zero",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "1",
              "maxcovarcount": "1",
              "rscape_covary_count": "0",
              "rnazprob": "0.95",
              "sci": "0.58",
              "consensus_mfe": "-14.10",
              "alifoldzscore": "-3.21",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )
    na_summary.write_text(
        textwrap.dedent(
            """\
            {
              "wbn": "rec_na",
              "nrseq": "2",
              "alilen": "12",
              "maxcovarval": "1",
              "maxcovarcount": "1",
              "rscape_covary_count": "NA",
              "rnazprob": "0.95",
              "sci": "0.58",
              "consensus_mfe": "-14.10",
              "alifoldzscore": "-3.21",
              "alifold_consstruc": "<<>>"
            }
            """
        ),
        encoding="utf-8",
    )

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "write-summary-outputs",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "RNAConSnake.log"),
            "--csv",
            str(tmp_path / "RNAConSnake.log.csv"),
            "--markdown",
            str(tmp_path / "RNAConSnake.md"),
            str(zero_summary),
            str(na_summary),
        ],
        check=True,
        env=subprocess_env(),
    )

    log_text = read_text(tmp_path / "RNAConSnake.log")
    csv_text = read_text(tmp_path / "RNAConSnake.log.csv")
    md_text = read_text(tmp_path / "RNAConSnake.md")

    assert "rec_zero" in log_text and "rscape 0" in log_text
    assert "rec_na" in log_text and "rscape NA" in log_text
    assert "rec_zero,2,12,1,1,0,0.95,0.58,-14.10,-3.21,<<>>" in csv_text
    assert "rec_na,2,12,1,1,NA,0.95,0.58,-14.10,-3.21,<<>>" in csv_text
    assert "| rec_zero | 2 | 12 | 1 | 1 | 0 | 0.95 | 0.58 | -14.10 | -3.21 | <<>> |" in md_text
    assert "| rec_na | 2 | 12 | 1 | 1 | NA | 0.95 | 0.58 | -14.10 | -3.21 | <<>> |" in md_text
    assert "<<>>" in log_text


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


_REFOLD_RECORD = re.compile(r"^> \S+\n[ACGUN]+\n[().]+ \( *-?\d+\.\d\d\)$", re.M)


def test_cli_workflow_smoke_test_with_fake_rnalalifold(tmp_path: Path) -> None:
    input_alignment = tmp_path / "my_input.stk"
    export_dir = tmp_path / "bundle"
    input_alignment.write_text(
        "# STOCKHOLM 1.0\n#=GF ID dummy_input\nseq1 ACGU\n//\n",
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
        "CLUSTAL W\n\nseq1 ACGU\nseq2 ACGU\n***\n",
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


def test_docs_reference_supported_input_formats() -> None:
    for doc in [
        Path("README.md"),
        Path("docs/usage.md"),
        Path("docs/pipeline_summary.md"),
    ]:
        text = read_text(doc)
        assert ".stk" in text
        if doc.name == "usage.md":
            assert "{stk,aln}" in text
        else:
            assert ".aln" in text


# --- null-model calibration arm ------------------------------------------


FAKE_RANDOMIZE_ALN = """\
#!/usr/bin/env perl
# Fake rnazRandomizeAln.pl for CI: permutes alignment columns, which preserves
# per-sequence gap counts exactly like the real column shuffler does.
use strict;
use warnings;

my @names;
my %seq;
while (my $line = <>) {
    chomp $line;
    next if $line =~ /^CLUSTAL/i;
    next if $line =~ /^\\s*$/;
    my ($name, $chunk) = split /\\s+/, $line, 2;
    next unless defined $chunk;
    push @names, $name unless exists $seq{$name};
    $seq{$name} .= $chunk;
}
die "no sequences\\n" unless @names;
my $len = length($seq{$names[0]});
my @order = (0 .. $len - 1);
for (my $i = $#order; $i > 0; $i--) {
    my $j = int(rand($i + 1));
    @order[$i, $j] = @order[$j, $i];
}
print "CLUSTAL W (1.81) multiple sequence alignment\\n\\n\\n";
for my $name (@names) {
    my $out = join '', map { substr($seq{$name}, $_, 1) } @order;
    print "$name $out\\n";
}
"""


def write_fake_randomize_aln(bin_dir: Path) -> Path:
    script = bin_dir / "rnazRandomizeAln.pl"
    script.write_text(FAKE_RANDOMIZE_ALN, encoding="utf-8")
    script.chmod(0o755)
    return script


TOY_ALIGNMENT = "\n".join(
    [
        "# STOCKHOLM 1.0",
        "#=GF ID toy",
        "seqA GGCUAGCUAGCUAACGUAGCUAGCUAGGCAUCGAUCGAUCG",
        "seqB GGCUAGCUAG---ACGUAGCUAGCUAGGCAUCGAUCGAUCG",
        "seqC GGCUAGCUAGCUAACGUAGCUCGCUAGGCAUCGAUCG---G",
        "//",
        "",
    ]
)


def test_null_settings_default_to_disabled_and_legacy_paths() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    disabled = NullSettings.from_config({})
    assert disabled.method == "none"
    assert disabled.enabled is False
    assert disabled.arms() == []

    enabled = NullSettings.from_config({"null": {"method": "sissiz", "replicates": 3}})
    assert enabled.enabled is True
    assert enabled.arms() == ["real", "null_000", "null_001", "null_002"]

    # replicates: 0 must be as inert as method: none
    assert NullSettings.from_config({"null": {"method": "sissiz", "replicates": 0}}).enabled is False


def test_null_settings_accepts_unquoted_yaml_null_key() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    # An unquoted "null:" key in YAML parses as the null scalar, so the section
    # arrives under None. It must not be silently ignored.
    settings = NullSettings.from_config({None: {"method": "sissiz", "replicates": 2}})
    assert settings.method == "sissiz"
    assert settings.replicates == 2


def test_null_settings_rejects_unknown_method() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    with pytest.raises(ValueError, match="Unknown null.method"):
        NullSettings.from_config({"null": {"method": "dinucleotide-magic"}})


def test_arm_seed_is_deterministic_arm_specific_and_none_for_real() -> None:
    from rnaconsnake.workflow_helpers import arm_seed

    assert arm_seed("real", 20261101) is None
    assert arm_seed("null_000", 20261101) == arm_seed("null_000", 20261101)
    assert arm_seed("null_000", 20261101) != arm_seed("null_001", 20261101)
    assert arm_seed("null_000", 20261101) != arm_seed("null_000", 20261102)


def test_candidate_paths_prefixed_by_arm_directory() -> None:
    plain = CandidatePaths(wlen=100, file="RC_100_0001_aln_1_12")
    armed = CandidatePaths(wlen=100, file="RC_100_0001_aln_1_12", arm_prefix="arms/null_000/")

    assert plain.rnaz_json == "generated_files/rnaz/len_100/RC_100_0001_aln_1_12.rnaz.json"
    assert armed.rnaz_json == "arms/null_000/" + plain.rnaz_json
    assert armed.aln_pdf == "arms/null_000/" + plain.aln_pdf
    assert armed.split == "arms/null_000/" + plain.split


def test_stockholm_parser_concatenates_interleaved_blocks(tmp_path: Path) -> None:
    from rnaconsnake.tools.stockholm_utils import parse_stockholm_records

    path = tmp_path / "interleaved.stk"
    path.write_text(
        "\n".join(
            [
                "# STOCKHOLM 1.0",
                "seqA ACGU",
                "seqB ACGA",
                "#=GC SS_cons <<>>",
                "",
                "seqA GGCC",
                "seqB GGCU",
                "#=GC SS_cons <<>>",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    record = parse_stockholm_records(path)[0]
    assert record.seqs["seqA"] == "ACGUGGCC"
    assert record.seqs["seqB"] == "ACGAGGCU"
    assert record.gc_lines == ["#=GC SS_cons <<>><<>>"]


def test_clustal_roundtrip_and_composition_diagnostics(tmp_path: Path) -> None:
    from rnaconsnake.tools.alignment_io import (
        alignment_diagnostics,
        parse_clustal_blocks,
        read_stockholm_alignment,
        write_clustal,
    )

    source = tmp_path / "toy.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    alignment = read_stockholm_alignment(source)
    assert alignment.order == ["seqA", "seqB", "seqC"]
    assert alignment.length == 41

    clustal = tmp_path / "toy.aln"
    write_clustal(alignment, clustal, width=20)
    blocks = parse_clustal_blocks(read_text(clustal))
    assert len(blocks) == 1
    assert blocks[0].seqs == alignment.seqs

    diagnostics = alignment_diagnostics(alignment)
    assert diagnostics["n_seq"] == 3
    assert diagnostics["length"] == 41
    assert 0.0 < diagnostics["mean_pairwise_identity"] <= 1.0
    assert abs(sum(diagnostics["base_composition"].values()) - 1.0) < 1e-6


def test_null_model_make_arm_copies_real_alignment_verbatim(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import make_arm_alignment

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    output = tmp_path / "arms" / "real" / "input_alignment.stk"
    make_arm_alignment("real", source, output)

    # A copy, never a symlink: downstream tools rewrite alignments in place.
    assert not output.is_symlink()
    assert read_text(output) == TOY_ALIGNMENT


def test_null_model_pool_roundtrip_with_fake_backend(tmp_path: Path) -> None:
    from rnaconsnake.tools.alignment_io import read_stockholm_alignment
    from rnaconsnake.tools.null_model import make_arm_alignment, simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover - perl is present on CI
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    metadata = simulate_pool(
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        method="rnazRandomizeAln",
        replicates=2,
        seed=20261101,
        sissiz_command=["SISSIz"],
        randomize_command=[str(randomize)],
        workdir=tmp_path / "work",
    )
    assert metadata["seeded"] is True
    assert metadata["determinism"] == "seeded"
    assert len(metadata["replicate_diagnostics"]) == 2

    original = read_stockholm_alignment(source)
    for arm in ["null_000", "null_001"]:
        out = tmp_path / arm / "input_alignment.stk"
        make_arm_alignment(arm, source, out, pool_path=tmp_path / "pool.stk")
        replicate = read_stockholm_alignment(out)
        assert replicate.order == original.order
        assert replicate.length == original.length
        # Column shuffling preserves per-sequence gap counts exactly.
        for name in original.order:
            assert replicate.seqs[name].count("-") == original.seqs[name].count("-")

    first = read_text(tmp_path / "null_000" / "input_alignment.stk")
    second = read_text(tmp_path / "null_001" / "input_alignment.stk")
    assert first != second


def test_null_model_pool_is_reproducible_for_the_same_seed(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)
    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    def run(tag: str, seed: int) -> str:
        simulate_pool(
            source_path=source,
            output_path=tmp_path / f"pool_{tag}.stk",
            metadata_path=tmp_path / f"pool_{tag}.json",
            method="rnazRandomizeAln",
            replicates=2,
            seed=seed,
            sissiz_command=["SISSIz"],
            randomize_command=[str(randomize)],
            workdir=tmp_path / f"work_{tag}",
        )
        return read_text(tmp_path / f"pool_{tag}.stk")

    assert run("a", 20261101) == run("b", 20261101)
    assert run("a", 20261101) != run("c", 20261102)


def test_null_model_rejects_single_sequence_alignment(tmp_path: Path) -> None:
    from rnaconsnake.tools.null_model import NullModelError, simulate_pool

    source = tmp_path / "input.stk"
    source.write_text("# STOCKHOLM 1.0\nseqA ACGU\n//\n", encoding="utf-8")
    with pytest.raises(NullModelError, match="at least 2 sequences"):
        simulate_pool(
            source_path=source,
            output_path=tmp_path / "pool.stk",
            metadata_path=tmp_path / "pool.json",
            method="rnazRandomizeAln",
            replicates=1,
            seed=1,
            sissiz_command=["SISSIz"],
            randomize_command=["rnazRandomizeAln.pl"],
        )


def _window_rows(spans: list[tuple[int, int]], structures: dict[str, str] | None = None):
    rows = []
    for start, end in spans:
        name = f"RC_100_0001_aln_{start}_{end}"
        rows.append(
            {
                "wbn": name,
                "rnazprob": "0.5",
                "alifoldzscore": "-1.0",
                "alifold_consstruc": (structures or {}).get(name, ""),
            }
        )
    return rows


def test_cluster_by_overlap_merges_any_shared_column() -> None:
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates
    from rnaconsnake.tools.loci import collapse_ratio

    rows = _window_rows([(1, 50), (40, 90), (91, 150), (200, 250)])
    clusters = cluster_candidates(
        candidates_from_records(rows), method="overlap", min_overlap=1, label="len100"
    )
    assert [(c.start, c.end, len(c.members)) for c in clusters] == [
        (1, 90, 2),
        (91, 150, 1),
        (200, 250, 1),
    ]
    assert [c.locus_id for c in clusters] == ["len100_0001", "len100_0002", "len100_0003"]
    assert collapse_ratio(4, 3) == pytest.approx(4 / 3)


def test_cluster_by_containment_keeps_adjacent_elements_apart() -> None:
    """Overlapping is not the same as nested; only fragments should collapse."""
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    # 712-799 with two nested fragments, plus 694-713 which merely overlaps it.
    rows = _window_rows([(712, 799), (719, 794), (745, 759), (694, 713)])
    clusters = cluster_candidates(candidates_from_records(rows), method="containment", label="len100")
    spans = [(c.start, c.end, len(c.members)) for c in clusters]
    assert spans == [(694, 713, 1), (712, 799, 3)]

    # Single-linkage overlap clustering would chain all four into one span.
    chained = cluster_candidates(
        candidates_from_records(rows), method="overlap", min_overlap=1, label="len100"
    )
    assert [(c.start, c.end, len(c.members)) for c in chained] == [(694, 799, 4)]


def test_cluster_none_leaves_every_window_alone() -> None:
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(1, 100), (11, 90), (21, 80)])
    clusters = cluster_candidates(candidates_from_records(rows), method="none", label="len100")
    assert len(clusters) == 3


def test_cluster_rejects_unknown_method() -> None:
    from rnaconsnake.tools.dereplicate import cluster_candidates

    with pytest.raises(ValueError, match="Unknown de-replication method"):
        cluster_candidates([], method="telepathy")


def test_base_pairs_map_to_absolute_alignment_columns() -> None:
    from rnaconsnake.tools.dereplicate import base_pairs

    assert base_pairs("((..))", offset=0) == {(0, 5), (1, 4)}
    assert base_pairs("((..))", offset=100) == {(100, 105), (101, 104)}
    # Unbalanced input must not raise; RNAalifold output is well formed but a
    # truncated record must not abort a run.
    assert base_pairs("((..", offset=0) == set()
    assert base_pairs("", offset=0) == set()


def test_substructure_clustering_uses_base_pair_subsets() -> None:
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    structures = {
        "RC_100_0001_aln_1_20": "((((....))))........",
        "RC_100_0001_aln_1_12": "((((....))))",
        "RC_100_0001_aln_5_16": "(((......)))",
    }
    rows = _window_rows([(1, 20), (1, 12), (5, 16)], structures)
    clusters = cluster_candidates(
        candidates_from_records(rows), method="substructure", pair_threshold=0.9, label="len100"
    )
    grouped = {c.locus_id: sorted(m.name for m in c.members) for c in clusters}
    # 1-12 shares every pair with 1-20; 5-16 folds differently and stays apart.
    assert any(
        set(members) == {"RC_100_0001_aln_1_12", "RC_100_0001_aln_1_20"} for members in grouped.values()
    )
    assert any(members == ["RC_100_0001_aln_5_16"] for members in grouped.values())


def test_representative_is_the_best_scoring_member() -> None:
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(1, 100), (11, 90)])
    rows[0]["rnazprob"] = "0.20"
    rows[1]["rnazprob"] = "0.97"
    clusters = cluster_candidates(candidates_from_records(rows), method="containment")
    assert len(clusters) == 1
    # The nested window scores better, so it represents the locus.
    assert clusters[0].representative.name == "RC_100_0001_aln_11_90"


def test_parse_window_name_rejects_unparseable_candidates() -> None:
    from rnaconsnake.tools.loci import parse_window_name

    with pytest.raises(ValueError, match="Cannot parse alignment coordinates"):
        parse_window_name("not_a_candidate_name")


def test_empirical_fdr_and_qvalue_envelope_are_monotone() -> None:
    from rnaconsnake.tools.calibration import empirical_fdr, qvalue_envelope

    real = [0.95, 0.90, 0.80, 0.60]
    nulls = [[0.60, 0.55], [0.62, 0.50]]
    fdr = empirical_fdr(real, nulls)
    envelope = qvalue_envelope(fdr)

    assert set(fdr) == {0.60, 0.80, 0.90, 0.95}
    # Nothing in the null arms reaches 0.80, so the high-score end is clean.
    assert fdr[0.95] == 0.0
    assert fdr[0.80] == 0.0
    assert fdr[0.60] == pytest.approx(1.0 / 4)

    values = [envelope[threshold] for threshold in sorted(envelope)]
    assert values == sorted(values, reverse=True) or all(
        values[i] >= values[i + 1] for i in range(len(values) - 1)
    )
    assert all(0.0 <= value <= 1.0 for value in values)


def _write_summary_csv(path: Path, rows: list[dict[str, str]]) -> str:
    from rnaconsnake.workflow_helpers import SUMMARY_FIELDS

    path.parent.mkdir(parents=True, exist_ok=True)
    import csv as _csv

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})
    return str(path)


def _candidate_rows(prob: float, alifoldz: float, rscape: str, count: int = 4) -> list[dict[str, str]]:
    rows = []
    for index in range(count):
        start = 1 + index * 500
        rows.append(
            {
                "wbn": f"RC_100_0001_aln_{start}_{start + 99}",
                "rnazprob": f"{prob:.4f}",
                "alifoldzscore": f"{alifoldz:.4f}",
                "rscape_covary_count": rscape,
                "nrseq": "6",
                "alilen": "100",
            }
        )
    return rows


def test_calibration_summary_records_every_clustering_parameter(tmp_path: Path) -> None:
    """Clustering decides how many loci each arm reports, and the q-values are
    counts over loci. A summary missing those parameters cannot be reproduced
    from -- and the export manifest copies this block verbatim."""
    from dataclasses import fields

    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
    }
    for index in range(2):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", _candidate_rows(0.20, 0.5, "0"))}

    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "rnazRandomizeAln", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )

    recorded = summary["thresholds"]
    for field in fields(Thresholds):
        assert field.name in recorded, f"{field.name} is not recorded in summary.json"
    assert recorded["max_container_width"] == 120
    assert recorded["container_min_coverage"] == pytest.approx(0.8)
    assert recorded["representative_rule"] == "widest"


def test_calibrate_writes_funnel_qvalues_and_summary(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
    }
    for index in range(3):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", _candidate_rows(0.20, 0.5, "0"))}

    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "rnazRandomizeAln", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )

    assert summary["counting_unit"] == "merged_loci"
    assert summary["fdr_conditional_on_stage_one"] is True
    assert summary["cascade_fdr"] == 0.0
    assert summary["counts"]["real_loci"] == 4
    assert summary["q_resolution"] == pytest.approx(1 / 3)
    assert summary["warnings"] == []

    funnel = read_text(tmp_path / "calibration" / "funnel.tsv")
    assert "# fdr_conditional_on_stage_one\ttrue" in funnel
    assert "counts are on merged loci" in funnel
    # One row per filter stage per arm, plus the null_mean aggregate.
    for arm in ["real", "null_000", "null_001", "null_002", "null_mean"]:
        for stage in ["windows", "loci", "rnaz", "alifoldz", "rscape", "cascade"]:
            assert f"{arm}\t100\t{stage}\t" in funnel

    qvalues = read_text(tmp_path / "calibration" / "qvalues.tsv")
    assert qvalues.count("\tyes\t") == 4
    assert "q_cascade" in qvalues

    dists = read_text(tmp_path / "calibration" / "score_dists.tsv")
    assert "null_000\t100\t" in dists
    assert "real\t" not in dists.split("\n", 1)[1]


def test_calibrate_is_reproducible_for_identical_inputs(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "2"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.30, -0.5, "0"))},
        "null_001": {100: _write_summary_csv(tmp_path / "n1.csv", _candidate_rows(0.95, -3.0, "1"))},
    }
    thresholds = Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2)
    for tag in ["a", "b"]:
        calibrate(
            arm_inputs=arm_inputs,
            thresholds=thresholds,
            null_metadata={"method": "sissiz", "seed": 20261101, "warnings": []},
            output_dir=tmp_path / tag,
            two_stage=True,
        )
    assert read_text(tmp_path / "a" / "qvalues.tsv") == read_text(tmp_path / "b" / "qvalues.tsv")
    assert read_text(tmp_path / "a" / "funnel.tsv") == read_text(tmp_path / "b" / "funnel.tsv")


def test_calibrate_warns_when_collapse_ratio_diverges_between_arms(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    # Real arm: 4 windows collapsing into 1 locus. Null arm: 4 windows, 4 loci.
    real_rows = [
        {
            "wbn": f"RC_100_0001_aln_{1 + index * 10}_{100 - index * 10}",
            "rnazprob": "0.95",
            "alifoldzscore": "-3.0",
            "rscape_covary_count": "1",
            "nrseq": "6",
            "alilen": "100",
        }
        for index in range(4)
    ]
    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", real_rows)},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=False,
    )
    assert summary["warnings"]
    assert "collapse ratio" in summary["warnings"][0]
    assert "# WARNING\t" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_calibrate_rejects_stage_one_looser_than_reported_threshold(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "1"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.1, 1.0, "0"))},
    }
    with pytest.raises(ValueError, match="stage1_rnaz_prob"):
        calibrate(
            arm_inputs=arm_inputs,
            thresholds=Thresholds(0.9, -2.0, 1, 0.95, 1, 0.2),
            null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
            output_dir=tmp_path / "calibration",
            two_stage=True,
        )


def test_calibrate_treats_stage_one_skips_as_failing_alifoldz(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    skipped = _candidate_rows(0.97, -4.0, "1")
    for row in skipped:
        row["alifoldzscore"] = "NA"
    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", skipped)},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.1, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    # A missing AlifoldZ score must never be read as a good (very negative) one.
    assert summary["counts"]["real_cascade_survivors"] == 0


BENCHMARK_TRUTH = "\n".join(
    [
        "element_id\telement_class\talignment\tstart\tend\tnotes",
        "xrRNA1\txrRNA\tflavivirus_3utr\t100\t199\tprimary target",
        "xrRNA2\txrRNA\tflavivirus_3utr\t900\t999\tsecond copy",
        "",
    ]
)


def test_benchmark_scores_recovery_against_calibrated_loci(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, read_qvalues, read_truth, write_recovery

    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(BENCHMARK_TRUTH, encoding="utf-8")
    qvalues_path = tmp_path / "qvalues.tsv"
    qvalues_path.write_text(
        "\n".join(
            [
                "# fdr_conditional_on_stage_one\ttrue",
                "locus_id\twlen\tstart\tend\trnazprob\tq_rnaz\tq_alifoldz\tq_cascade\tcascade_pass",
                "len100_0001\t100\t90\t210\t0.97\t0\t0.01\t0\tyes",
                "len100_0002\t100\t400\t500\t0.5\t0.4\t0.5\tNA\tno",
                "",
            ]
        ),
        encoding="utf-8",
    )

    results = evaluate(
        read_truth(truth_path, "flavivirus_3utr"),
        read_qvalues(qvalues_path),
        min_overlap_fraction=0.5,
        allow_uncurated=False,
    )
    by_id = {row["element_id"]: row for row in results}
    assert by_id["xrRNA1"]["recovered"] == "yes"
    assert by_id["xrRNA1"]["best_locus"] == "len100_0001"
    assert by_id["xrRNA1"]["q_cascade"] == "0"
    assert by_id["xrRNA2"]["recovered"] == "no"

    out = tmp_path / "recovery.tsv"
    write_recovery(results, out)
    text = read_text(out)
    assert "# recovered\t1" in text
    assert "# curated\t2" in text


def test_benchmark_refuses_uncurated_truth_file(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth

    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "xrRNA1\txrRNA\tflavivirus_3utr\tTBD\tTBD\tnot curated yet\n",
        encoding="utf-8",
    )
    truth = read_truth(truth_path)
    with pytest.raises(BenchmarkError, match="uncurated coordinates"):
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)

    results = evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=True)
    assert results[0]["recovered"] == "uncurated"


def test_shipped_flavivirus_truth_file_is_schema_valid() -> None:
    from rnaconsnake.tools.benchmark import TRUTH_COLUMNS, read_truth

    truth = read_truth(Path("resources/benchmark/flavivirus_elements.tsv"))
    assert {"xrRNA", "DB", "sHP", "3SL"} <= {element.element_class for element in truth}
    assert TRUTH_COLUMNS[0] == "element_id"
    # Coordinates are alignment specific and intentionally not shipped curated.
    assert all(not element.curated for element in truth)


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


def test_check_deps_requires_the_viennarna_python_module(monkeypatch) -> None:
    """The bindings are as much a runtime dependency as the binaries, and the
    preflight is where a missing one has to surface."""
    monkeypatch.setattr(cli, "viennarna_bindings_version", lambda: None)
    assert cli.check_dependencies() == 2


def _fake_viennarna_tools(tmp_path: Path, version: str) -> dict[str, str]:
    """Stand-ins for every ViennaRNA binary the check looks at, so the machine's
    own install does not decide the outcome."""
    tools = {}
    for key in cli.VIENNARNA_BINARIES:
        script = tmp_path / key
        script.write_text(f"#!/usr/bin/env python3\nprint('{key} {version}')\n", encoding="utf-8")
        script.chmod(0o755)
        tools[key] = str(script)
    return tools


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
    assert "ViennaRNA version mismatch" in captured.err
    assert "same ViennaRNA release" in captured.err


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


def test_snakefile_defines_every_rule_exactly_once() -> None:
    import re as _re

    text = read_text(Path("snakefile"))
    names = _re.findall(r"^\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", text, _re.M)
    duplicates = sorted({name for name in names if names.count(name) > 1})
    assert not duplicates, f"duplicated rule definitions: {duplicates}"

    # No parallel arm-specific rules: the arms must traverse the same DAG.
    assert not [name for name in names if name.endswith(("_null", "_real"))]
    assert "make_arm_alignment" in names
    assert "calibrate" in names


def test_snakefile_has_no_duplicated_rule_bodies() -> None:
    import re as _re

    text = read_text(Path("snakefile"))
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.M)
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
    assert null_0 != real_arm and null_1 != real_arm and null_0 != null_1

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
    assert _re.search(r"^make_arm_alignment\s+3$", output, _re.M), output
    assert _re.search(r"^calibrate\s+1$", output, _re.M), output
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


def test_calibrate_drops_rscape_from_cascade_when_it_never_ran(tmp_path: Path) -> None:
    """do_rscape: false makes every count NA; that must not zero the headline."""
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "NA"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "NA"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    assert summary["rscape_evaluated"] is False
    assert summary["cascade_filters"] == ["rnaz", "alifoldz"]
    assert summary["counts"]["real_cascade_survivors"] == 4
    assert summary["cascade_fdr"] == 0.0
    assert any("R-scape produced no covariation counts" in w for w in summary["warnings"])
    assert "# rscape_in_cascade\tfalse" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_calibrate_keeps_rscape_in_cascade_when_it_ran(tmp_path: Path) -> None:
    from rnaconsnake.tools.calibration import Thresholds, calibrate

    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", _candidate_rows(0.97, -4.0, "0"))},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.10, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    assert summary["rscape_evaluated"] is True
    assert summary["cascade_filters"] == ["rnaz", "alifoldz", "rscape"]
    # Zero covarying pairs is a real, reported result, not a missing value.
    assert summary["counts"]["real_cascade_survivors"] == 0
    assert "# rscape_in_cascade\ttrue" in read_text(tmp_path / "calibration" / "funnel.tsv")


def test_split_passthrough_config_merges_user_entries() -> None:
    """Snakemake's --config is nargs="+", so two occurrences lose the first."""
    rest, entries = cli.split_passthrough_config(
        ["--configfile", "c.yaml", "--config", "a=1", "b=2", "--cores", "4"]
    )
    assert rest == ["--configfile", "c.yaml", "--cores", "4"]
    assert entries == ["a=1", "b=2"]

    assert cli.split_passthrough_config(["--dry-run"]) == (["--dry-run"], [])
    assert cli.split_passthrough_config(["-C", "x=1"]) == ([], ["x=1"])


# --- lower-case alignment regression --------------------------------------


LOWERCASE_CLUSTAL = "\n".join(
    [
        "CLUSTAL 2.1 multiple sequence alignment",
        "",
        "NMV_NC_032088.1 aggcacagaacgccg",
        "KOKV_NC_009029. cggcacagaacgccg",
        "ZIKV_NC_012532. aggcacagatcgccg",
        "                *********:*****",
        "",
    ]
)


def test_clean_clustal_uppercases_sequences_for_refold(tmp_path: Path) -> None:
    """The cleaned Clustal is a recorded artifact, and it stays upper-case.

    It was refold.pl that forced this -- it matches [A-Z-] rows and read a
    lower-case alignment as nothing at all. The Python refold that replaced it
    reads any case, but the artifact contract has not changed.
    """
    source = tmp_path / "in.aln"
    source.write_text(LOWERCASE_CLUSTAL, encoding="utf-8")
    backup = tmp_path / "in.aln~"
    output = tmp_path / "cleaned.aln"

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "clean-clustal",
            "--input",
            str(source),
            "--backup",
            str(backup),
            "--output",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    cleaned = read_text(output)
    assert "NMV_NC_032088.1 AGGCACAGAACGCCG" in cleaned
    assert "aggcacagaacgccg" not in cleaned
    # Names keep their case, the header is untouched, conservation rows go away.
    assert "CLUSTAL 2.1 multiple sequence alignment" in cleaned
    assert "*" not in cleaned
    # The backup is the untouched original.
    assert read_text(backup) == LOWERCASE_CLUSTAL

    for line in cleaned.splitlines():
        if line and not line.startswith("CLUSTAL"):
            assert not any(char.islower() for char in line.split(" ", 1)[-1])


def test_clean_clustal_is_a_noop_for_uppercase_alignments(tmp_path: Path) -> None:
    upper = LOWERCASE_CLUSTAL.replace("aggcacagaacgccg", "AGGCACAGAACGCCG")
    upper = upper.replace("cggcacagaacgccg", "CGGCACAGAACGCCG")
    upper = upper.replace("aggcacagatcgccg", "AGGCACAGATCGCCG")
    source = tmp_path / "in.aln"
    source.write_text(upper, encoding="utf-8")
    output = tmp_path / "cleaned.aln"

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "clean-clustal",
            "--input",
            str(source),
            "--backup",
            str(tmp_path / "b.aln~"),
            "--output",
            str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    expected = "\n".join(line for line in upper.splitlines() if "*" not in line).rstrip() + "\n"
    assert read_text(output) == expected


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


# --- alifoldz sentinel / case-normalisation regression --------------------


ALIFOLDZ_EMPTY_REPORT = """\
###################################################################
# alifoldz.pl
#
#           Input: 0 sequences of 0 columns
#   Sample Number: 100
###################################################################

  From      To    Strand    Native MFE    Mean MFE     STDV      Z
 ------------------------------------------------------------------

9999
"""

ALIFOLDZ_GOOD_REPORT = """\
###################################################################
#           Input: 20 sequences of 52 columns
###################################################################

  From      To    Strand    Native MFE    Mean MFE     STDV      Z
 ------------------------------------------------------------------
     1      52       +         -4.64       -1.77       1.61    -1.8
-1.8
"""


def _run_extract_alifoldz(tmp_path: Path, report: str):
    source = tmp_path / "a.alifoldz.txt"
    source.write_text(report, encoding="utf-8")
    output = tmp_path / "a.alifoldz.json"
    return (
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.legacy_postprocess",
                "extract-alifoldz",
                "--input",
                str(source),
                "--output",
                str(output),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=subprocess_env(),
        ),
        output,
    )


def test_extract_alifoldz_rejects_empty_input_report(tmp_path: Path) -> None:
    """alifoldz exits 0 having parsed nothing; 9999 is its 'no window' sentinel."""
    result, output = _run_extract_alifoldz(tmp_path, ALIFOLDZ_EMPTY_REPORT)
    assert result.returncode != 0
    assert "0 sequences" in result.stderr
    assert "lower-case" in result.stderr
    assert not output.exists()


def test_extract_alifoldz_rejects_the_9999_sentinel(tmp_path: Path) -> None:
    report = ALIFOLDZ_EMPTY_REPORT.replace(
        "Input: 0 sequences of 0 columns", "Input: 6 sequences of 40 columns"
    )
    result, output = _run_extract_alifoldz(tmp_path, report)
    assert result.returncode != 0
    assert "9999" in result.stderr
    assert not output.exists()


def test_extract_alifoldz_accepts_a_real_score(tmp_path: Path) -> None:
    result, output = _run_extract_alifoldz(tmp_path, ALIFOLDZ_GOOD_REPORT)
    assert result.returncode == 0, result.stderr
    assert json.loads(read_text(output)) == {"alifoldzscore": "-1.8"}


def test_strip_aln_uppercases_sequences_for_the_perl_toolchain(tmp_path: Path) -> None:
    """alifoldz.pl matches [A-Z...] with no /i flag, and is still Perl."""
    source = tmp_path / "lower.stk"
    source.write_text(
        "# STOCKHOLM 1.0\n"
        "#=GF ID lower\n"
        "Seq_a tcagaccacttag--tgcc\n"
        "Seq_b tcagatccgaaa---ggcc\n"
        "#=GC SS_cons <<<<.......>>>>....\n"
        "//\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.tools.strip_aln", "-a", str(source), "-f", "S", "--nosingle"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert "Seq_a TCAGACCACTTAG--TGCC" in result.stdout
    assert "Seq_b TCAGATCCGAAA---GGCC" in result.stdout
    # Names and annotation rows keep their original form.
    assert "Seq_a" in result.stdout and "SEQ_A" not in result.stdout
    assert "#=GC SS_cons <<<<.......>>>>...." in result.stdout


def test_strip_aln_leaves_uppercase_alignments_untouched(tmp_path: Path) -> None:
    source = tmp_path / "upper.stk"
    body = "# STOCKHOLM 1.0\n#=GF ID upper\nSeq_a TCAGACCACTTAG--TGCC\nSeq_b TCAGATCCGAAA---GGCC\n//\n"
    source.write_text(body, encoding="utf-8")
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.tools.strip_aln", "-a", str(source), "-f", "S", "--nosingle"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert "Seq_a TCAGACCACTTAG--TGCC" in result.stdout
    assert "Seq_b TCAGATCCGAAA---GGCC" in result.stdout


def test_strip_aln_dedups_case_variant_duplicates(tmp_path: Path) -> None:
    """Normalising before dedup means case variants collapse, as they should."""
    source = tmp_path / "mixed.stk"
    source.write_text(
        "# STOCKHOLM 1.0\nSeq_a acgu\nSeq_b ACGU\nSeq_c ACGA\n//\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [PYTHON, "-m", "rnaconsnake.tools.strip_aln", "-a", str(source), "-f", "S"],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert "Seq_a ACGU" in result.stdout
    assert "Seq_b" not in result.stdout
    assert "Seq_c ACGA" in result.stdout


def test_non_redundant_table_reports_one_row_per_locus(tmp_path: Path) -> None:
    """The de-replicated table keeps the best window and records what it absorbed."""
    from rnaconsnake.tools.dereplicate import (
        NR_COLUMNS,
        candidates_from_records,
        cluster_candidates,
        write_non_redundant_csv,
    )

    rows = [
        {
            "wbn": "RC_100_0001_aln_712_799",
            "rnazprob": "0.999984",
            "alifoldzscore": "-2.3",
            "alifold_consstruc": "",
            "nrseq": "20",
            "alilen": "88",
        },
        {
            "wbn": "RC_100_0001_aln_739_778",
            "rnazprob": "0.969731",
            "alifoldzscore": "-0.9",
            "alifold_consstruc": "",
            "nrseq": "20",
            "alilen": "40",
        },
        {
            "wbn": "RC_100_0001_aln_694_713",
            "rnazprob": "0.888226",
            "alifoldzscore": "-3.8",
            "alifold_consstruc": "",
            "nrseq": "20",
            "alilen": "20",
        },
    ]
    clusters = cluster_candidates(candidates_from_records(rows), method="containment", label="len100")
    output = tmp_path / "RNAConSnake.nr.csv"
    write_non_redundant_csv(clusters, output)

    import csv as _csv

    with open(output, encoding="utf-8") as handle:
        written = list(_csv.DictReader(handle))
    assert [row["locus_id"] for row in written] == ["len100_0002", "len100_0001"]
    assert list(written[0]) == NR_COLUMNS

    by_locus = {row["locus_id"]: row for row in written}
    absorbed = by_locus["len100_0002"]
    assert absorbed["wbn"] == "RC_100_0001_aln_712_799"
    assert absorbed["n_windows"] == "2"
    # Nothing is silently dropped: the fragment is named in the members column.
    assert "RC_100_0001_aln_739_778" in absorbed["members"]

    # The merely-overlapping neighbour keeps its own row.
    adjacent = by_locus["len100_0001"]
    assert adjacent["wbn"] == "RC_100_0001_aln_694_713"
    assert adjacent["n_windows"] == "1"


def test_calibration_counts_the_same_loci_the_ranked_table_reports(tmp_path: Path) -> None:
    """q-values must describe the objects that actually get reported."""
    from rnaconsnake.tools.calibration import Thresholds, calibrate
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    # One real element reported at its extent plus three nested fragments.
    spans = [(1, 100), (11, 90), (21, 80), (31, 70)]
    rows = [
        {
            "wbn": f"RC_100_0001_aln_{s}_{e}",
            "rnazprob": "0.97",
            "alifoldzscore": "-4.0",
            "rscape_covary_count": "2",
            "nrseq": "6",
            "alilen": str(e - s + 1),
            "alifold_consstruc": "",
        }
        for s, e in spans
    ]
    arm_inputs = {
        "real": {100: _write_summary_csv(tmp_path / "real.csv", rows)},
        "null_000": {100: _write_summary_csv(tmp_path / "n0.csv", _candidate_rows(0.1, 1.0, "0"))},
    }
    summary = calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2, "containment", 0.9),
        null_metadata={"method": "sissiz", "seed": 1, "warnings": []},
        output_dir=tmp_path / "calibration",
        two_stage=True,
    )
    assert summary["dereplication_method"] == "containment"
    assert summary["counts"]["real_windows"] == 4
    # Four windows, one element: the FDR is computed on one locus, not four.
    assert summary["counts"]["real_loci"] == 1

    clusters = cluster_candidates(candidates_from_records(rows), method="containment")
    assert len(clusters) == summary["counts"]["real_loci"]


def test_config_declares_a_dereplication_method() -> None:
    import yaml as _yaml

    from rnaconsnake.tools.dereplicate import METHODS

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    assert payload["dereplicate"]["method"] in METHODS


def test_configurable_rules_declare_their_settings_as_params() -> None:
    """Config read from globals inside `run:` is invisible to Snakemake.

    Without a `params:` declaration the params rerun-trigger cannot see a
    changed setting, so re-running with a different method silently reuses the
    previous output.
    """
    import re as _re

    text = read_text(Path("snakefile"))
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.M)
    by_name = {}
    for block in blocks:
        match = _re.match(r"\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", block)
        if match:
            by_name[match.group(1)] = block

    for rule in ["dereplicate_summary", "calibrate", "benchmark_recovery"]:
        assert rule in by_name, f"{rule} is missing from the snakefile"
        assert _re.search(r"^\s+params:", by_name[rule], _re.M), (
            f"rule {rule} reads configurable settings but declares no params:, "
            "so changing them will not trigger a re-run"
        )

    # ...and the run body must consume params, not the globals directly,
    # otherwise the declaration is decorative.
    body = by_name["dereplicate_summary"].split("\n    run:", 1)[1]
    assert "params.method" in body
    assert "DEREPLICATE_METHOD" not in body
    assert "PAIR_CONTAINMENT" not in body


def test_dereplication_doc_matches_the_implementation() -> None:
    """The spec in docs/dereplication.md must not drift from the code."""
    from rnaconsnake.tools.dereplicate import METHODS, NR_COLUMNS
    from rnaconsnake.workflow_helpers import CalibrationSettings

    doc = read_text(Path("docs/dereplication.md"))
    defaults = CalibrationSettings.from_config({})

    # Every method is documented, and the documented default is the real one.
    for method in METHODS:
        assert f"`{method}`" in doc, f"method {method} is undocumented"
    assert f"`{defaults.dereplicate_method}` is the default" in doc

    # Documented column order matches what is written.
    assert "locus_id  locus_start  locus_end  n_windows  members" in doc
    assert NR_COLUMNS[:5] == ["locus_id", "locus_start", "locus_end", "n_windows", "members"]

    # Documented default thresholds match the shipped ones.
    assert f"default {defaults.pair_containment}" in doc
    assert f"default {defaults.locus_min_overlap}" in doc

    # The doc claims the step is unconditional, so the target must be among the
    # default targets rule all requests.
    snakefile = read_text(Path("snakefile"))
    targets = snakefile.split("def pipeline_targets():")[1].split("\nrule ")[0]
    assert "RNAConSnake.nr.csv" in targets
    assert "pipeline_targets()" in snakefile.split("rule all:")[1].split("\nrule ")[0]

    # And that the standalone invocation it prints actually exists.
    assert "python -m rnaconsnake.tools.dereplicate" in doc
    for flag in ["--input", "--output", "--method", "--label"]:
        assert flag in doc


def test_dereplication_doc_is_linked_from_the_user_docs() -> None:
    for path, needle in [
        (Path("README.md"), "docs/dereplication.md"),
        (Path("docs/usage.md"), "dereplication.md"),
        (Path("docs/pipeline_summary.md"), "dereplication.md"),
    ]:
        assert needle in read_text(path), f"{path} does not link the de-replication spec"


def _render_markdown(tmp_path: Path, nr_rows, full_rows, method="containment"):
    import csv as _csv

    from rnaconsnake.tools.dereplicate import NR_COLUMNS
    from rnaconsnake.workflow_helpers import SUMMARY_FIELDS

    nr = tmp_path / "RNAConSnake.nr.csv"
    with open(nr, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=NR_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in nr_rows:
            writer.writerow({field: row.get(field, "") for field in NR_COLUMNS})

    full = tmp_path / "RNAConSnake.log.csv"
    with open(full, "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in full_rows:
            writer.writerow({field: row.get(field, "") for field in SUMMARY_FIELDS})

    out = tmp_path / "RNAConSnake.md"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "render-markdown",
            "--label",
            "len_200",
            "--nr",
            str(nr),
            "--full",
            str(full),
            "--output",
            str(out),
            "--method",
            method,
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    return read_text(out)


def test_markdown_puts_the_non_redundant_block_before_the_full_table(tmp_path: Path) -> None:
    nr_rows = [
        {
            "locus_id": "len200_0006",
            "locus_start": "717",
            "locus_end": "803",
            "n_windows": "2",
            "members": "RC_200_0001_aln_717_803;RC_200_0001_aln_723_799",
            "wbn": "RC_200_0001_aln_717_803",
            "rnazprob": "0.999987",
            "alifoldzscore": "-2.2",
        }
    ]
    full_rows = [
        {"wbn": "RC_200_0001_aln_717_803", "rnazprob": "0.999987", "alifoldzscore": "-2.2"},
        {"wbn": "RC_200_0001_aln_723_799", "rnazprob": "0.999968", "alifoldzscore": "-1.8"},
    ]
    text = _render_markdown(tmp_path, nr_rows, full_rows)

    assert text.startswith("# RNAConSnake Summary: len_200")
    assert "- Loci (non-redundant): 1" in text
    assert "- Windows (all): 2" in text
    assert "- De-replication method: `containment`" in text

    nr_at = text.index("## Non-redundant candidates")
    full_at = text.index("## All windows")
    assert nr_at < full_at, "the non-redundant block must come first"

    # Both blocks are present in one file; there is no separate .md.
    assert "len200_0006" in text[nr_at:full_at]
    assert "RC_200_0001_aln_723_799" in text[full_at:]
    assert "including the 1 collapsed above" in text

    # members is pushed to the last column: widest and least often read.
    header = [line for line in text.splitlines() if line.startswith("| locus_id")][0]
    assert header.rstrip().endswith("members |")


def test_markdown_escapes_pipes_so_tables_do_not_break(tmp_path: Path) -> None:
    nr_rows = [{"locus_id": "len200_0001", "wbn": "w1", "alifold_consstruc": "((|))"}]
    text = _render_markdown(tmp_path, nr_rows, [{"wbn": "w1", "alifold_consstruc": "((|))"}])
    assert r"((\|))" in text


def test_summary_logs_no_longer_owns_the_markdown() -> None:
    """The Markdown report is rendered downstream of de-replication.

    summary_logs cannot render it: the leading block needs the non-redundant
    table, which is produced from summary_logs' own CSV.
    """
    import re as _re

    text = read_text(Path("snakefile"))
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.M)
    by_name = {}
    for block in blocks:
        match = _re.match(r"\s*(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", block)
        if match:
            by_name[match.group(1)] = block

    assert "RNAConSnake.md" not in by_name["summary_logs"]
    renderer = by_name["render_summary_markdown"]
    assert "RNAConSnake.md" in renderer
    assert "RNAConSnake.nr.csv" in renderer
    assert "RNAConSnake.log.csv" in renderer


def test_write_summary_outputs_markdown_is_optional(tmp_path: Path) -> None:
    summary = tmp_path / "one.summary.json"
    summary.write_text('{"wbn": "rec1", "rnazprob": "0.95"}\n', encoding="utf-8")
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.legacy_postprocess",
            "write-summary-outputs",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "RNAConSnake.log"),
            "--csv",
            str(tmp_path / "RNAConSnake.log.csv"),
            str(summary),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert (tmp_path / "RNAConSnake.log").is_file()
    assert not (tmp_path / "RNAConSnake.md").exists()


# --- export bundle: de-replication and calibration annotations -------------


def test_export_annotations_resolve_the_real_arm(tmp_path: Path) -> None:
    """With the null arm on, pipeline outputs live under arms/real/."""
    from rnaconsnake.export_annotations import analysis_root

    plain = tmp_path / "plain"
    (plain / "generated_files" / "summary" / "len_100").mkdir(parents=True)
    assert analysis_root(plain) == plain

    armed = tmp_path / "armed"
    (armed / "arms" / "real" / "generated_files" / "summary" / "len_100").mkdir(parents=True)
    assert analysis_root(armed) == armed / "arms" / "real"


def _write_nr_table(run_dir: Path, wlen: int, rows: list[dict[str, str]], method="containment"):
    import csv as _csv

    from rnaconsnake.tools.dereplicate import NR_COLUMNS

    target = run_dir / "generated_files" / "summary" / f"len_{wlen}"
    target.mkdir(parents=True, exist_ok=True)
    with open(target / "RNAConSnake.nr.csv", "w", encoding="utf-8", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=NR_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in NR_COLUMNS})
    (target / "RNAConSnake.nr.json").write_text(
        json.dumps({"method": method, "n_loci": len(rows)}) + "\n", encoding="utf-8"
    )


def test_export_annotations_map_windows_to_loci(tmp_path: Path) -> None:
    from rnaconsnake.export_annotations import collect

    _write_nr_table(
        tmp_path,
        100,
        [
            {
                "locus_id": "len100_0001",
                "locus_start": "712",
                "locus_end": "799",
                "n_windows": "2",
                "members": "RC_100_0001_aln_712_799;RC_100_0001_aln_739_778",
                "wbn": "RC_100_0001_aln_712_799",
            }
        ],
    )
    annotations = collect(tmp_path)
    assert annotations.has_loci is True
    assert annotations.dereplication_method == "containment"
    assert annotations.locus_count == 1
    assert annotations.representative_count == 1

    representative = annotations.locus_for("RC_100_0001_aln_712_799")
    fragment = annotations.locus_for("RC_100_0001_aln_739_778")
    assert representative is fragment, "both windows belong to the same locus"
    assert representative.is_representative("RC_100_0001_aln_712_799") is True
    assert representative.is_representative("RC_100_0001_aln_739_778") is False

    # No calibration present: the bundle must still be describable.
    assert annotations.has_calibration is False
    assert annotations.qvalues_for("RC_100_0001_aln_712_799") is None


def test_export_annotations_read_calibration_qvalues(tmp_path: Path) -> None:
    from rnaconsnake.export_annotations import collect, manifest_block

    _write_nr_table(
        tmp_path,
        100,
        [
            {
                "locus_id": "len100_0001",
                "locus_start": "1",
                "locus_end": "100",
                "n_windows": "1",
                "members": "RC_100_0001_aln_1_100",
                "wbn": "RC_100_0001_aln_1_100",
            }
        ],
    )
    calibration = tmp_path / "results" / "calibration"
    calibration.mkdir(parents=True)
    (calibration / "qvalues.tsv").write_text(
        "# fdr_conditional_on_stage_one\ttrue\n"
        "locus_id\tq_rnaz\tq_alifoldz\tq_cascade\tcascade_pass\n"
        "len100_0001\t0\t0.25\t0.1\tyes\n",
        encoding="utf-8",
    )
    (calibration / "summary.json").write_text(
        json.dumps(
            {
                "null": {"method": "sissiz", "replicates": 100, "seed": 1},
                "two_stage": True,
                "fdr_conditional_on_stage_one": True,
                "fdr_conditionality_note": "conditional",
                "cascade_fdr": 0.1,
                "cascade_filters": ["rnaz", "alifoldz"],
                "counting_unit": "merged_loci",
                "warnings": [],
            }
        ),
        encoding="utf-8",
    )

    annotations = collect(tmp_path)
    assert annotations.has_calibration is True
    q = annotations.qvalues_for("RC_100_0001_aln_1_100")
    assert q is not None and q.q_cascade == "0.1" and q.cascade_pass == "yes"

    block = manifest_block(annotations)
    assert block["calibration"]["applied"] is True
    assert block["calibration"]["null_method"] == "sissiz"
    assert block["calibration"]["replicates"] == 100
    # The conditionality of the FDR must be legible without reading any table.
    assert block["calibration"]["fdr_conditional_on_stage_one"] is True


def test_export_bundle_columns_cover_loci_and_qvalues() -> None:
    from rnaconsnake.export_bundle import CANDIDATE_COLUMNS, FEATURE_COLUMNS

    for column in [
        "locus_id",
        "locus_start",
        "locus_end",
        "locus_window_count",
        "is_representative",
        "redundant_to",
        "q_rnaz",
        "q_alifoldz",
        "q_cascade",
        "cascade_pass",
    ]:
        assert column in CANDIDATE_COLUMNS
    assert "locus_count" in FEATURE_COLUMNS
    assert "representative_count" in FEATURE_COLUMNS


def test_locus_columns_are_blank_without_dereplication() -> None:
    """A bundle from an uncalibrated run stays readable by the same consumer."""
    from rnaconsnake.export_annotations import RunAnnotations
    from rnaconsnake.export_bundle import locus_columns

    columns = locus_columns(RunAnnotations(), "RC_100_0001_aln_1_100")
    assert set(columns) == {
        "locus_id",
        "locus_start",
        "locus_end",
        "locus_window_count",
        "is_representative",
        "redundant_to",
        "q_rnaz",
        "q_alifoldz",
        "q_cascade",
        "cascade_pass",
    }
    assert all(value == "" for value in columns.values())


# --- scope boundary: computation, not presentation --------------------------
#
# RNAConSnake produces alignments, figures, metrics and structured export
# bundles. Rendering or browsing those results is out of scope. No templates,
# vendored web assets, rendering code, or web-framework dependency belongs in
# this repository; the export bundle is the interface downstream consumers read.

WEB_FORBIDDEN_PATTERNS = [
    r"\bjinja2?\b",
    r"\bdatatables\b",
    r"\btabler\b",
    r"<!doctype html",
    r"<html[ >]",
    r"render_.*_pages",
    r"site_builder",
]

WEB_FORBIDDEN_PATHS = (".html", ".jinja", ".j2", ".css", ".js")


def _tracked_files() -> list[Path]:
    result = subprocess.run(["git", "ls-files"], check=True, capture_output=True, text=True, cwd=str(ROOT))
    return [Path(line) for line in result.stdout.splitlines() if line]


def test_no_web_assets_are_tracked() -> None:
    offenders = [path for path in _tracked_files() if path.suffix.lower() in WEB_FORBIDDEN_PATHS]
    assert not offenders, "web assets are out of scope for RNAConSnake: " + ", ".join(
        str(path) for path in offenders
    )


def test_no_rendering_code_or_dependency_is_tracked() -> None:
    import re as _re

    pattern = _re.compile("|".join(WEB_FORBIDDEN_PATTERNS), _re.IGNORECASE)
    offenders: list[str] = []
    for path in _tracked_files():
        full = ROOT / path
        if not full.is_file():
            continue
        # This test names the forbidden things, so it cannot police itself.
        if full.resolve() == Path(__file__).resolve():
            continue
        try:
            text = full.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for match in pattern.finditer(text):
            offenders.append(f"{path}: {match.group(0)!r}")
    assert not offenders, "presentation-layer code leaked into RNAConSnake:\n" + "\n".join(offenders)


def test_packaging_declares_no_web_dependencies() -> None:
    text = read_text(Path("pyproject.toml")).lower()
    for forbidden in ["jinja", "flask", "fastapi", "markdown", "starlette"]:
        assert forbidden not in text, f"{forbidden} is a presentation concern, not a pipeline one"


def test_workflow_emits_no_html() -> None:
    """AGENTS.md: no browser-style HTML reporting in this project."""
    snakefile = read_text(Path("snakefile"))
    assert ".html" not in snakefile
    for path in _tracked_files():
        assert not str(path).endswith(".html")


def test_sissiz_backend_is_reported_as_unseeded() -> None:
    """SISSIz draws its seed from the clock, so null.seed cannot reproduce it.

    Two invocations inside the same second agree, which makes this easy to
    mis-measure; the metadata must not claim reproducibility it cannot deliver.
    """
    import inspect

    from rnaconsnake.tools import null_model

    source = inspect.getsource(null_model._simulate_sissiz)
    assert '"determinism": "unseeded"' in source
    assert '"seeded": False' in source
    assert "tool_is_deterministic" not in source


def test_unseeded_backend_warns_that_the_seed_does_not_reproduce_the_pool(
    tmp_path: Path,
) -> None:
    from rnaconsnake.tools.null_model import simulate_pool

    if shutil.which("perl") is None:  # pragma: no cover
        pytest.skip("perl is required for the rnazRandomizeAln backend")

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    randomize = write_fake_randomize_aln(bin_dir)
    # Make the fake report itself as a non-Perl command so it takes the
    # unseeded path, exercising the warning.
    plain = bin_dir / "randomize_plain"
    plain.write_text(randomize.read_text(encoding="utf-8"), encoding="utf-8")
    plain.chmod(0o755)

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    metadata = simulate_pool(
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        method="rnazRandomizeAln",
        replicates=1,
        seed=1,
        sissiz_command=["SISSIz"],
        randomize_command=[str(plain)],
        workdir=tmp_path / "work",
    )
    assert metadata["seeded"] is False
    assert any("not seedable" in warning for warning in metadata["warnings"])


def test_adopt_pool_reuses_a_pinned_pool_and_validates_it(tmp_path: Path) -> None:
    """Pinning is how a SISSIz calibration is made reproducible."""
    from rnaconsnake.tools.null_model import NullModelError, adopt_pool

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")

    pool = tmp_path / "pinned.stk"
    pool.write_text(
        "\n".join(
            [
                "# STOCKHOLM 1.0",
                "#=GF ID null_000",
                "seqA GGCUAGCUAGCUAACGUAGCUAGCUAGGCAUCGAUCGAUCG",
                "seqB GGCUAGCUAG---ACGUAGCUAGCUAGGCAUCGAUCGAUCG",
                "seqC GGCUAGCUAGCUAACGUAGCUCGCUAGGCAUCGAUCG---G",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )

    metadata = adopt_pool(
        pool_path=pool,
        source_path=source,
        output_path=tmp_path / "pool.stk",
        metadata_path=tmp_path / "pool.json",
        replicates=1,
        method="sissiz",
    )
    assert metadata["determinism"] == "pinned"
    assert metadata["pinned_from"] == str(pool)
    assert read_text(tmp_path / "pool.stk") == read_text(pool)

    with pytest.raises(NullModelError, match="need 2"):
        adopt_pool(
            pool_path=pool,
            source_path=source,
            output_path=tmp_path / "pool2.stk",
            metadata_path=tmp_path / "pool2.json",
            replicates=2,
            method="sissiz",
        )


def test_adopt_pool_rejects_a_pool_from_a_different_alignment(tmp_path: Path) -> None:
    """A pinned pool from another alignment would silently miscalibrate."""
    from rnaconsnake.tools.null_model import NullModelError, adopt_pool

    source = tmp_path / "input.stk"
    source.write_text(TOY_ALIGNMENT, encoding="utf-8")
    wrong = tmp_path / "wrong.stk"
    wrong.write_text(
        "# STOCKHOLM 1.0\n#=GF ID null_000\nseqA ACGU\nseqB ACGA\nseqC ACGC\n//\n",
        encoding="utf-8",
    )
    with pytest.raises(NullModelError, match="length"):
        adopt_pool(
            pool_path=wrong,
            source_path=source,
            output_path=tmp_path / "pool.stk",
            metadata_path=tmp_path / "pool.json",
            replicates=1,
            method="sissiz",
        )


def test_perl_seeded_command_wraps_scripts_and_falls_back(tmp_path: Path) -> None:
    from rnaconsnake.workflow_helpers import perl_seed_env, perl_seeded_command

    script = tmp_path / "helper.pl"
    script.write_text("print 1;\n", encoding="utf-8")
    command, seeded = perl_seeded_command([str(script)], 42, ["-f"])
    assert seeded is True
    assert command[0] == "perl" and command[3] == "42" and command[-1] == "-f"

    # Non-Perl commands run unchanged, and say so rather than pretending.
    fallback, seeded = perl_seeded_command(["/bin/echo"], 42, ["hi"])
    assert seeded is False
    assert fallback == ["/bin/echo", "hi"]

    with pytest.raises(ValueError):
        perl_seeded_command([], 1)

    # Hash order must be pinned too; srand alone is not enough.
    env = perl_seed_env({})
    assert env["PERL_HASH_SEED"] == "0"
    assert env["PERL_PERTURB_KEYS"] == "0"


def test_derived_seed_is_deterministic_and_per_candidate() -> None:
    from rnaconsnake.workflow_helpers import derived_seed

    assert derived_seed(7, "cand_a") == derived_seed(7, "cand_a")
    assert derived_seed(7, "cand_a") != derived_seed(7, "cand_b")
    assert derived_seed(7, "cand_a") != derived_seed(8, "cand_a")


def test_snakefile_seeds_alifoldz_per_candidate() -> None:
    """alifoldz.pl shuffles internally and has no seed option.

    Its z-scores feed the q-values, so an unseeded AlifoldZ makes the whole
    calibration irreproducible no matter how the null pool was generated.
    """
    text = read_text(Path("snakefile"))
    assert 'ALIFOLDZ_SEED = config.get("alifoldz_seed")' in text
    assert "perl_seeded_command(" in text
    assert "derived_seed(int(ALIFOLDZ_SEED), wildcards.file)" in text


def test_config_seeds_alifoldz_by_default() -> None:
    import yaml as _yaml

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    assert payload.get("alifoldz_seed") is not None, "a calibrated run must be reproducible by default"


def test_null_settings_accepts_a_pinned_pool() -> None:
    from rnaconsnake.workflow_helpers import NullSettings

    settings = NullSettings.from_config({"null": {"method": "sissiz", "pool_file": "/tmp/pool.stk"}})
    assert settings.pool_file == "/tmp/pool.stk"
    # Snakemake stringifies nested --config values, so "None" must not be a path.
    for sentinel in ["None", "", "null", "~"]:
        assert (
            NullSettings.from_config({"null": {"method": "sissiz", "pool_file": sentinel}}).pool_file is None
        )


# --- sensitivity: long windows must not swallow distinct elements ----------


def test_wide_window_does_not_absorb_small_distinct_elements() -> None:
    """A 200 nt RNALalifold window over three 60 nt hairpins is a scan
    artifact, not their parent. Absorbing them reports one locus for three."""
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(1, 200), (10, 70), (90, 150), (160, 200)])
    merged = cluster_candidates(candidates_from_records(rows), method="containment", max_container_width=0)
    assert len(merged) == 1, "without the guard everything collapses"

    guarded = cluster_candidates(candidates_from_records(rows), method="containment", max_container_width=120)
    assert len(guarded) == 4


def test_wide_windows_still_absorb_near_duplicates_of_themselves() -> None:
    """Otherwise every long overlapping window becomes its own locus."""
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(1, 200), (4, 198), (6, 195)])
    clusters = cluster_candidates(
        candidates_from_records(rows),
        method="containment",
        max_container_width=120,
        min_coverage=0.8,
    )
    assert len(clusters) == 1


def test_narrow_container_still_absorbs_its_fragments() -> None:
    """The 3'SL case: an 88 nt element keeps absorbing its own sub-hairpins."""
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(712, 799), (719, 794), (739, 778), (745, 759)])
    clusters = cluster_candidates(
        candidates_from_records(rows), method="containment", max_container_width=120
    )
    assert len(clusters) == 1
    assert len(clusters[0].members) == 4


def test_alignment_report_separates_gappy_from_low_identity(tmp_path: Path) -> None:
    """Gappiness prevents detection; low identity only costs power. Conflating
    them would flag regions where strong elements are in fact still found."""
    from rnaconsnake.tools.alignment_io import Alignment
    from rnaconsnake.tools.alignment_report import profile, summarize

    # First half: only one of four sequences has residues. Second half: full
    # occupancy but highly divergent.
    seqs = {
        "a": "ACGUACGUAC" + "ACGUACGUAC",
        "b": "----------" + "GCAUGCAUGC",
        "c": "----------" + "UAUAUAUAUA",
        "d": "----------" + "CGCGCGCGCG",
    }
    alignment = Alignment(order=list(seqs), seqs=seqs)
    regions = profile(alignment, window=10, min_occupancy=0.5, min_identity=0.55)
    assert regions[0].screenable is False
    assert any("occupancy" in reason for reason in regions[0].reasons)
    assert regions[1].screenable is True
    assert regions[1].reduced_power is True

    payload = summarize(regions, alignment)
    assert payload["blind_columns"] == 10
    assert payload["reduced_power_columns"] == 10
    assert any("not evidence of absence" in w for w in payload["warnings"])


def test_alignment_report_is_a_default_workflow_target() -> None:
    text = read_text(Path("snakefile"))
    targets = text.split("def pipeline_targets():")[1].split("\nrule ")[0]
    assert "alignment_screenability.tsv" in targets


# --- benchmark truth scaffolding -------------------------------------------


def test_structural_domains_descend_through_long_range_pairs() -> None:
    """Long-range helices bracket several elements rather than defining one."""
    from rnaconsnake.tools.benchmark_scaffold import structural_domains

    # One long-range pair enclosing two element-scale hairpins.
    left = "(((((....)))))"  # 14 nt, 5 bp
    right = "((((......))))"  # 14 nt, 5 bp
    structure = "(" + left + "....." + right + ")"

    wide = structural_domains(structure, max_width=200, min_width=5, min_pairs=3)
    assert len(wide) == 1, "with a generous width the whole thing is one domain"

    narrow = structural_domains(structure, max_width=20, min_width=5, min_pairs=3)
    assert [(d.start, d.end) for d in narrow] == [(2, 15), (21, 34)]
    assert [d.n_pairs for d in narrow] == [5, 4]


def test_structural_domains_skip_trivial_helices() -> None:
    from rnaconsnake.tools.benchmark_scaffold import structural_domains

    assert structural_domains("((...))", max_width=50, min_width=5, min_pairs=3) == []
    assert structural_domains(".........", max_width=50) == []


def test_read_ss_cons_requires_a_reference_structure(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark_scaffold import read_ss_cons

    bare = tmp_path / "bare.stk"
    bare.write_text("# STOCKHOLM 1.0\nseqA ACGU\nseqB ACGA\n//\n", encoding="utf-8")
    with pytest.raises(ValueError, match="no '#=GC SS_cons' line"):
        read_ss_cons(bare)

    annotated = tmp_path / "annotated.stk"
    annotated.write_text("# STOCKHOLM 1.0\nseqA ACGU\nseqB ACGA\n#=GC SS_cons <<>>\n//\n", encoding="utf-8")
    structure, n_seq = read_ss_cons(annotated)
    assert structure == "<<>>" and n_seq == 2


def test_scaffold_marks_missed_domains_and_unexplained_loci() -> None:
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows

    domains = [Domain(start=10, end=60, n_pairs=12), Domain(start=200, end=260, n_pairs=14)]
    loci = [
        {"locus_id": "L1", "locus_start": "5", "locus_end": "70", "wbn": "w1", "rnazprob": "0.99"},
        {"locus_id": "L2", "locus_start": "400", "locus_end": "460", "wbn": "w2", "rnazprob": "0.97"},
        # Low-scoring and unmatched: not worth a curator's attention.
        {"locus_id": "L3", "locus_start": "600", "locus_end": "660", "wbn": "w3", "rnazprob": "0.1"},
    ]
    rows = scaffold_rows(domains, loci, "demo")
    by_id = {row["element_id"]: row for row in rows}

    assert "NOT REPORTED" not in by_id["TBD_01"]["notes"]
    assert "NOT REPORTED" in by_id["TBD_02"]["notes"], "a domain with no locus must be flagged"
    assert "TBD_extra_L2" in by_id, "a strong locus outside the reference needs a decision"
    assert "TBD_extra_L3" not in by_id
    # Every row is coordinate-complete but label-free.
    assert all(row["element_id"].startswith("TBD") for row in rows)
    assert all(row["start"] and row["end"] for row in rows)


def test_scaffold_does_not_double_report_sibling_loci() -> None:
    """The same element is reported once per window length."""
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows

    domains = [Domain(start=10, end=60, n_pairs=12)]
    loci = [
        {"locus_id": "len100_0001", "locus_start": "8", "locus_end": "62", "wbn": "a", "rnazprob": "0.99"},
        {"locus_id": "len200_0001", "locus_start": "5", "locus_end": "70", "wbn": "b", "rnazprob": "0.98"},
    ]
    rows = scaffold_rows(domains, loci, "demo")
    assert len(rows) == 1
    assert not any("extra" in row["element_id"] for row in rows)


def test_scaffold_output_is_readable_by_the_benchmark(tmp_path: Path) -> None:
    """The scaffold must be schema-compatible, and refused until curated."""
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth
    from rnaconsnake.tools.benchmark_scaffold import Domain, scaffold_rows, write_scaffold

    path = tmp_path / "scaffold.tsv"
    write_scaffold(scaffold_rows([Domain(10, 60, 12)], [], "demo"), path, "demo")

    truth = read_truth(path)
    assert len(truth) == 1
    # Coordinates are filled in, so it is the *labels* that block it, and the
    # benchmark still refuses it until a curator has been through.
    assert truth[0].curated is True
    assert truth[0].element_id.startswith("TBD")
    with pytest.raises(BenchmarkError, match="placeholder"):
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)


def test_shipped_jevg_truth_file_is_fully_curated() -> None:
    """The JEV-group truth file is curated repository data, not a scaffold."""
    from rnaconsnake.tools.benchmark import BenchmarkError, evaluate, read_truth

    truth = read_truth(Path("resources/benchmark/jevg_3utr_elements.tsv"), "jevg_3utr")
    assert len(truth) == 9
    assert all(element.curated for element in truth)

    names = {element.element_id for element in truth}
    assert {"xrRNA1", "xrRNA2", "DB1", "DB2", "sHP", "3SL"} <= names
    # HP1 and HP2 were not previously annotated; they are part of the finding.
    assert {"HP1", "HP2"} <= names

    # Spans are ordered, non-overlapping and within the alignment.
    spans = sorted((e.start, e.end) for e in truth)
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:], strict=False))
    assert spans[0][0] >= 1 and spans[-1][1] <= 711

    # No placeholders survive, so the benchmark accepts it.
    try:
        evaluate(truth, [], min_overlap_fraction=0.5, allow_uncurated=False)
    except BenchmarkError as error:  # pragma: no cover - guards a regression
        raise AssertionError(f"curated truth file was refused: {error}") from error


def test_flavivirus_truth_file_remains_an_uncurated_stub() -> None:
    """The broad MBFV file is still a stub; it must not be mistaken for curated."""
    from rnaconsnake.tools.benchmark import read_truth

    truth = read_truth(Path("resources/benchmark/flavivirus_elements.tsv"))
    assert all(not element.curated for element in truth)


# --- cascade threshold sweep ------------------------------------------------


def _sweep_arm_inputs(tmp_path: Path, real_rows, null_rows_per_arm):
    arm_inputs = {"real": {100: _write_summary_csv(tmp_path / "real.csv", real_rows)}}
    for index, rows in enumerate(null_rows_per_arm):
        arm = f"null_{index:03d}"
        arm_inputs[arm] = {100: _write_summary_csv(tmp_path / f"{arm}.csv", rows)}
    return arm_inputs


def _sweep_base():
    from rnaconsnake.tools.calibration import Thresholds

    return Thresholds(0.9, -2.0, 1, 0.5, 1, 0.2, "containment", 0.9, 120, 0.8)


def test_threshold_sweep_reports_survivors_and_fdr(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import sweep

    strong = _candidate_rows(0.97, -4.0, "NA")
    weak = _candidate_rows(0.20, 0.5, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, strong, [weak, weak]),
        rnaz_grid=[0.5, 0.9],
        alifoldz_grid=[-1.0, -3.0],
        base=_sweep_base(),
    )
    assert len(points) == 4
    by_key = {(p.rnaz_prob, p.alifoldz): p for p in points}

    # The real arm clears every combination; the null arm clears none.
    for point in points:
        assert point.real_survivors == 4
        assert point.null_mean == 0.0
        assert point.fdr == 0.0
    assert by_key[(0.9, -3.0)].null_sd == 0.0


def test_threshold_sweep_tightening_reduces_null_survivors(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.97, -4.0, "NA")
    # Null loci that clear a loose threshold but not a strict one.
    null = _candidate_rows(0.92, -1.2, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [null, null]),
        rnaz_grid=[0.9],
        alifoldz_grid=[-1.0, -3.0],
        base=_sweep_base(),
    )
    loose = [p for p in points if p.alifoldz == -1.0][0]
    strict = [p for p in points if p.alifoldz == -3.0][0]
    assert loose.null_mean > strict.null_mean
    assert loose.fdr > strict.fdr


def test_threshold_sweep_counts_reference_recovery(tmp_path: Path) -> None:
    """Recovery is sensitivity only, and must be reported alongside FDR."""
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.97, -4.0, "NA")  # spans 1-100, 501-600, 1001-1100, 1501-1600
    null = _candidate_rows(0.10, 1.0, "NA")
    spans = [(10, 90), (510, 590), (9000, 9100)]  # third is never covered
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [null]),
        rnaz_grid=[0.9],
        alifoldz_grid=[-2.0],
        base=_sweep_base(),
        reference_spans=spans,
    )
    assert points[0].recovered == 2
    assert points[0].reference_total == 3


def test_threshold_sweep_stage_one_gate_follows_the_reported_threshold(tmp_path: Path) -> None:
    """stage1 must never exceed the reported RNAz cutoff, or the sweep would
    select candidates whose AlifoldZ was never computed."""
    from rnaconsnake.tools.threshold_sweep import sweep

    real = _candidate_rows(0.6, -4.0, "NA")
    points = sweep(
        _sweep_arm_inputs(tmp_path, real, [_candidate_rows(0.1, 1.0, "NA")]),
        rnaz_grid=[0.2],
        alifoldz_grid=[-2.0],
        base=_sweep_base(),  # stage1 default is 0.5, above the 0.2 grid point
    )
    # No exception, and the loose threshold is honoured.
    assert points[0].real_survivors == 4


def test_write_sweep_marks_undefined_fdr(tmp_path: Path) -> None:
    from rnaconsnake.tools.threshold_sweep import SweepPoint, write_sweep

    output = tmp_path / "sweep.tsv"
    write_sweep(
        [SweepPoint(0.9, -2.0, 0, 0.0, 0.0, None, None, None)],
        output,
    )
    text = read_text(output)
    # No real survivors means no FDR, and that must not read as zero.
    assert "\tNA\tNA\tNA" in text
    assert "sensitivity only" in text


# --- benchmark null baseline ------------------------------------------------
#
# Reported loci cover most of an alignment, so an overlap test can be satisfied
# by chance. A recovery count without a null baseline is uninterpretable.


def _locus(locus_id, start, end, prob="0.99"):
    return {
        "locus_id": locus_id,
        "start": str(start),
        "end": str(end),
        "rnazprob": prob,
        "alifoldzscore": "-3.0",
        "q_cascade": "",
        "cascade_pass": "",
    }


def _two_elements(tmp_path: Path):
    from rnaconsnake.tools.benchmark import read_truth

    path = tmp_path / "truth.tsv"
    path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "e1\txrRNA\taln\t100\t199\tfirst\n"
        "e2\tDB\taln\t400\t499\tsecond\n",
        encoding="utf-8",
    )
    return read_truth(path, "aln")


def test_null_baseline_summarises_chance_recovery(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import null_baseline

    truth = _two_elements(tmp_path)
    arms = [
        [_locus("n1", 90, 210)],  # covers e1 only
        [_locus("n1", 90, 210), _locus("n2", 390, 510)],  # covers both
        [],  # covers neither
    ]
    baseline = null_baseline(truth, arms, min_overlap_fraction=0.5)
    assert baseline["arms"] == 3
    assert baseline["min"] == 0
    assert baseline["max"] == 2
    assert baseline["mean"] == pytest.approx(1.0)
    assert null_baseline(truth, [], 0.5) == {}


def test_recovery_report_warns_when_null_matches_the_real_arm(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    # Null arms recover more than the real arm did.
    write_recovery(results, output, {"arms": 5, "mean": 2.0, "min": 2, "max": 2})
    text = read_text(output)
    assert "# null_baseline_recovered\tmean 2.0" in text
    assert "recovery_margin_over_null\t-1.0" in text
    assert "not evidence of detection" in text


def test_recovery_report_warns_when_the_overlap_test_is_vacuous(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210), _locus("r2", 390, 510)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    write_recovery(results, output, {"arms": 10, "mean": 1.8, "min": 1, "max": 2})
    text = read_text(output)
    assert "recovery_margin_over_null\t+0.2" in text
    assert "close to vacuous" in text


def test_recovery_report_is_quiet_when_the_margin_is_real(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import evaluate, write_recovery

    truth = _two_elements(tmp_path)
    results = evaluate(truth, [_locus("r1", 90, 210), _locus("r2", 390, 510)], 0.5, allow_uncurated=True)
    output = tmp_path / "recovery.tsv"
    write_recovery(results, output, {"arms": 10, "mean": 0.4, "min": 0, "max": 1})
    text = read_text(output)
    assert "recovery_margin_over_null\t+1.6" in text
    assert "WARNING" not in text


def test_one_null_loci_group_is_one_arm(tmp_path: Path) -> None:
    """An arm's window lengths pool into a single baseline sample.

    Counting each window length as its own arm would dilute the mean with
    partial views of the same arm, and understate the baseline.
    """
    truth_path = tmp_path / "truth.tsv"
    truth_path.write_text(
        "element_id\telement_class\talignment\tstart\tend\tnotes\n"
        "e1\txrRNA\taln\t100\t199\tfirst\n"
        "e2\tDB\taln\t400\t499\tsecond\n",
        encoding="utf-8",
    )
    qvalues = tmp_path / "qvalues.tsv"
    qvalues.write_text(
        "locus_id\twlen\tstart\tend\trnazprob\talifoldzscore\nlen100_0001\t100\t90\t210\t0.99\t-3.0\n",
        encoding="utf-8",
    )
    header = "locus_id,locus_start,locus_end,n_windows,members,wbn,rnazprob,alifoldzscore\n"
    len100 = tmp_path / "len100.nr.csv"
    len100.write_text(header + "len100_0001,90,210,1,a,a,0.99,-3.0\n", encoding="utf-8")
    len200 = tmp_path / "len200.nr.csv"
    len200.write_text(header + "len200_0001,390,510,1,b,b,0.99,-3.0\n", encoding="utf-8")

    output = tmp_path / "recovery.tsv"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.benchmark",
            "--truth",
            str(truth_path),
            "--qvalues",
            str(qvalues),
            "--output",
            str(output),
            "--alignment",
            "aln",
            "--null-loci",
            str(len100),
            str(len200),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    text = read_text(output)
    # One arm, recovering both elements -- not two arms recovering one each.
    assert "over 1 null arms" in text
    assert "# null_baseline_recovered\tmean 2.0" in text


def test_benchmark_rule_passes_every_null_arm_as_one_group() -> None:
    """The recovery table is uninterpretable without the null baseline, so the
    rule that builds it has to supply the arms' locus tables itself."""
    text = read_text(Path("snakefile"))
    rule = text.split("rule benchmark_recovery:", 1)[1].split("\nrule ", 1)[0]
    assert "null_loci=" in rule
    assert 'cmd += ["--null-loci", *params.null_loci_by_arm[arm]]' in rule
    assert "for arm in sorted(params.null_loci_by_arm)" in rule


def test_reciprocal_overlap_exposes_an_oversized_locus(tmp_path: Path) -> None:
    """A locus far larger than the element scores 1.0 on overlap_fraction.

    That is how a screen reporting one huge locus can look perfectly sensitive.
    """
    from rnaconsnake.tools.benchmark import RECOVERY_COLUMNS, evaluate

    assert "reciprocal_overlap" in RECOVERY_COLUMNS
    truth = _two_elements(tmp_path)
    # One locus spanning the whole alignment contains both elements entirely.
    results = evaluate(truth, [_locus("huge", 1, 1000)], 0.5, allow_uncurated=True)
    for row in results:
        assert row["recovered"] == "yes"
        assert row["overlap_fraction"] == "1.0000"
        # ...but the element accounts for only a tenth of the locus.
        assert float(row["reciprocal_overlap"]) == pytest.approx(0.1, abs=1e-3)


def test_read_locus_table_accepts_the_non_redundant_csv(tmp_path: Path) -> None:
    from rnaconsnake.tools.benchmark import _read_locus_table

    path = tmp_path / "nr.csv"
    path.write_text(
        "locus_id,locus_start,locus_end,n_windows,members,wbn,rnazprob,alifoldzscore\n"
        "len100_0001,10,90,2,a;b,a,0.97,-3.1\n",
        encoding="utf-8",
    )
    rows = _read_locus_table(path)
    assert rows == [
        {"locus_id": "len100_0001", "start": "10", "end": "90", "rnazprob": "0.97", "alifoldzscore": "-3.1"}
    ]


def test_representative_rules_behave_as_documented() -> None:
    from rnaconsnake.tools.dereplicate import (
        REPRESENTATIVE_RULES,
        candidates_from_records,
        select_representative,
    )

    rows = _window_rows([(1, 100), (20, 60), (30, 45)])
    rows[0]["rnazprob"], rows[0]["alifoldzscore"] = "0.70", "-1.0"  # widest, weakest RNAz
    rows[1]["rnazprob"], rows[1]["alifoldzscore"] = "0.99", "-2.0"  # best RNAz
    rows[2]["rnazprob"], rows[2]["alifoldzscore"] = "0.80", "-4.0"  # best AlifoldZ
    members = candidates_from_records(rows)

    assert select_representative(members, "best_rnaz").width == 41
    assert select_representative(members, "widest").width == 100
    assert select_representative(members, "best_alifoldz").width == 16

    with pytest.raises(ValueError, match="Unknown representative rule"):
        select_representative(members, "vibes")
    assert set(REPRESENTATIVE_RULES) == {"best_rnaz", "widest", "best_alifoldz", "widest_of_top_half"}


def test_widest_representative_keeps_a_multi_lobed_element_whole() -> None:
    """The DB1 defect: a dumbbell represented by one lobe fails a cascade the
    whole element passes."""
    from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates

    rows = _window_rows([(425, 497), (431, 472)])
    rows[0]["rnazprob"], rows[0]["alifoldzscore"] = "0.992", "-2.4"  # whole element
    rows[1]["rnazprob"], rows[1]["alifoldzscore"] = "0.995", "-1.8"  # one lobe, scores better
    members = candidates_from_records(rows)

    by_rnaz = cluster_candidates(members, method="containment", representative_rule="best_rnaz")
    by_width = cluster_candidates(members, method="containment", representative_rule="widest")
    assert by_rnaz[0].representative.width == 42
    assert by_width[0].representative.width == 73
    # Only the wider representative clears alifoldz <= -2.0.
    assert float(by_rnaz[0].representative.record["alifoldzscore"]) > -2.0
    assert float(by_width[0].representative.record["alifoldzscore"]) <= -2.0


def test_configuration_is_locked_to_the_recorded_values() -> None:
    """The configuration was fixed on JEVG before DENVG was evaluated.

    Changing these values silently would forfeit the only out-of-sample number
    the study has; docs/configuration_lock.md records the reasoning.
    """
    import yaml as _yaml

    payload = _yaml.safe_load(read_text(Path("config.yaml")))
    assert payload["dereplicate"]["representative"] == "widest"
    assert payload["calibration"]["rnaz_prob_threshold"] == 0.8
    assert payload["calibration"]["alifoldz_threshold"] == -1.5

    lock = read_text(Path("docs/configuration_lock.md"))
    assert "representative: widest" in lock
    assert "rnaz_prob_threshold: 0.8" in lock
    assert "alifoldz_threshold: -1.5" in lock
    # The limitations must stay stated, not quietly dropped.
    assert "DB2 is not recovered" in lock
    assert "sensitive to a single element" in lock


def test_every_curated_truth_file_is_packaged() -> None:
    """`benchmark_truth` names a file in resources/benchmark by name, and an
    installed run resolves it inside the package. A truth file left out of the
    build hook exists in the repository and nowhere else."""
    import importlib.util
    import types

    # setup.py is read for its copy list alone, so setuptools is stubbed rather
    # than required: the test env does not build the package.
    setuptools = types.ModuleType("setuptools")
    setuptools.setup = lambda **kwargs: None
    command = types.ModuleType("setuptools.command")
    build_py = types.ModuleType("setuptools.command.build_py")
    build_py.build_py = type("build_py", (), {"run": lambda self: None})
    stubs = {
        "setuptools": setuptools,
        "setuptools.command": command,
        "setuptools.command.build_py": build_py,
    }
    saved = {name: sys.modules.get(name) for name in stubs}
    sys.modules.update(stubs)
    try:
        spec = importlib.util.spec_from_file_location("_rnaconsnake_setup", Path("setup.py"))
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        for name, previous in saved.items():
            if previous is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = previous

    packaged = {source.name for source in module.WORKFLOW_SOURCES}
    for truth in sorted(Path("resources/benchmark").glob("*.tsv")):
        assert truth.name in packaged, f"{truth.name} is not copied into the package"


# --- sensitivity envelope ---------------------------------------------------
#
# The null arm bounds false positives but says nothing about the detection
# floor, so "no candidates" on a divergent alignment is not evidence of absence
# until the floor is measured.


def _envelope_alignment():
    from rnaconsnake.tools.alignment_io import Alignment

    # Two near-identical sequences and two divergent ones, so subsets span a
    # range of mean pairwise identity.
    seqs = {
        "a": "ACGUACGUACGUACGUACGU",
        "b": "ACGUACGUACGUACGUACGA",
        "c": "AGGUAGGUAGGUAGGUAGGU",
        "d": "UCAUUCAUUCAUUCAUUCAU",
    }
    return Alignment(order=list(seqs), seqs=seqs)


def test_enumerate_subsets_is_exhaustive_when_cheap() -> None:
    from rnaconsnake.tools.sensitivity_envelope import enumerate_subsets

    alignment = _envelope_alignment()
    subsets = enumerate_subsets(alignment, n_seq=2, max_subsets=100)
    assert len(subsets) == 6  # C(4,2)
    # Sorted by identity, and the near-identical pair is the most similar.
    identities = [s.identity for s in subsets]
    assert identities == sorted(identities)
    assert set(subsets[-1].names) == {"a", "b"}


def test_enumerate_subsets_samples_deterministically_when_large() -> None:
    from rnaconsnake.tools.sensitivity_envelope import enumerate_subsets

    alignment = _envelope_alignment()
    first = enumerate_subsets(alignment, n_seq=2, max_subsets=3, seed=7)
    second = enumerate_subsets(alignment, n_seq=2, max_subsets=3, seed=7)
    assert len(first) == 3
    assert [s.names for s in first] == [s.names for s in second]
    assert [s.names for s in first] != [s.names for s in enumerate_subsets(alignment, 2, 3, seed=8)]


def test_enumerate_subsets_rejects_impossible_sizes() -> None:
    from rnaconsnake.tools.sensitivity_envelope import enumerate_subsets

    alignment = _envelope_alignment()
    assert enumerate_subsets(alignment, n_seq=1) == []
    assert enumerate_subsets(alignment, n_seq=99) == []


def test_pick_by_identity_spans_the_range_without_repeats() -> None:
    from rnaconsnake.tools.sensitivity_envelope import enumerate_subsets, pick_by_identity

    alignment = _envelope_alignment()
    subsets = enumerate_subsets(alignment, n_seq=2, max_subsets=100)
    chosen = pick_by_identity(subsets, [0.0, 0.5, 1.0])
    assert len(chosen) == 3
    assert len({c.names for c in chosen}) == 3, "must not return the same subset twice"
    # Asking for more targets than there are subsets stops rather than repeating.
    assert len(pick_by_identity(subsets, [0.5] * 20)) == len(subsets)


def test_plan_writes_subsets_that_preserve_columns(tmp_path: Path) -> None:
    """Subsetting drops sequences, never columns: element coordinates survive."""
    from rnaconsnake.tools.alignment_io import read_stockholm_alignment
    from rnaconsnake.tools.sensitivity_envelope import Subset, plan, write_subset_alignment

    alignment = _envelope_alignment()
    rows = plan(alignment, sizes=[2, 3], targets=[0.4, 0.8], max_subsets=100, seed=1)
    assert rows and all(row["label"].startswith("n") for row in rows)

    subset = Subset(names=tuple(rows[0]["names"].split(";")), identity=rows[0]["identity"])
    path = tmp_path / "subset.stk"
    write_subset_alignment(alignment, subset, path, "demo")
    written = read_stockholm_alignment(path)
    assert written.length == alignment.length
    assert written.n_seq == subset.n_seq


def test_shipped_denvg_truth_file_is_the_held_out_test_set() -> None:
    """DENVG was curated after the configuration was locked on JEVG.

    Four elements only: the scaffold domains that turned out not to be xrRNAs
    are excluded rather than labelled uncertainly, because an unbelieved row
    would inflate recovery.
    """
    from rnaconsnake.tools.benchmark import read_truth

    truth = read_truth(Path("resources/benchmark/denvg_3utr_elements.tsv"), "denvg_3utr")
    assert len(truth) == 4
    assert {e.element_id for e in truth} == {"DB1", "DB2", "sHP", "3SL"}
    assert all(e.curated for e in truth)

    spans = sorted((e.start, e.end) for e in truth)
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:], strict=False))
    assert spans[-1][1] <= 488

    text = read_text(Path("resources/benchmark/denvg_3utr_elements.tsv"))
    # The exclusions must stay documented, not silently dropped.
    assert "EXCLUDED" in text
    assert "not xrRNAs" in text
    assert "HELD-OUT" in text


def test_scaffold_flags_domains_that_need_attention(capsys) -> None:
    """Wide domains bracket several elements; narrow ones are easy to miss.

    Lowering max_width was tested and does not improve which elements are
    captured -- a wide domain already contains them -- so the scaffold warns
    rather than changing the decomposition.
    """
    from rnaconsnake.tools.benchmark_scaffold import DEFAULT_MIN_WIDTH, Domain, scaffold_rows

    # A confirmed 14 nt element was dropped by the previous min-width of 15.
    assert DEFAULT_MIN_WIDTH <= 14

    scaffold_rows([Domain(1, 131, 30), Domain(200, 214, 4)], [], "demo")
    # scaffold_rows itself is silent; the warnings live in main(). Assert the
    # thresholds those warnings use are the ones documented.
    import inspect

    from rnaconsnake.tools import benchmark_scaffold

    source = inspect.getsource(benchmark_scaffold.main)
    assert "exceed 100 nt" in source
    assert "under 20 nt" in source


def test_readme_documents_the_tools_conda_cannot_supply() -> None:
    """Three tools are on no package index; a user who does not know that is
    stuck before the first run."""
    readme = read_text(Path("README.md"))
    container = read_text(Path("container/README.md"))

    assert "https://github.com/mtw/SISSIz" in readme
    assert "https://github.com/mtw/SISSIz" in container
    for tool in ["SISSIz", "alifoldz.pl"]:
        assert tool in readme, f"{tool} is not mentioned in the README"
    # The --use-conda limitation is a stated submission blocker; it must not
    # silently disappear from the docs.
    assert "--use-conda" in readme

    # prepare-context.sh must name the sources in its failure path, since that
    # is where a user actually hits the problem.
    prepare = read_text(Path("container/prepare-context.sh"))
    assert "github.com/mtw/SISSIz" in prepare
    assert "RNAz source tarball" in prepare


def test_readme_covers_the_features_added_this_cycle() -> None:
    readme = read_text(Path("README.md"))
    for topic in [
        "null-model calibration",  # the calibration arm
        "De-replicating",  # locus de-replication
        "screenability",  # blind-region reporting
        "Positive control",  # benchmark
        "container",  # containerised toolchain
    ]:
        assert topic.lower() in readme.lower(), f"README does not mention {topic}"


def test_readme_links_resolve() -> None:
    """Links must point at files the repository actually ships.

    Checking only that a path exists locally is not enough: a file present in
    the working tree but untracked (or gitignored) resolves here and 404s for
    anyone who clones.
    """
    import re as _re

    readme = read_text(Path("README.md"))
    targets = [t for t in _re.findall(r"\[`[^`]+`\]\(([^)]+)\)", readme) if not t.startswith("http")]
    for target in targets:
        assert Path(target).exists(), f"README links to a missing path: {target}"

    # An unstaged file is merely not added yet; an *ignored* one will never be
    # in a clone, so the link is permanently broken for everyone else.
    ignored = subprocess.run(
        ["git", "check-ignore", *targets],
        capture_output=True,
        text=True,
        cwd=str(ROOT),
    ).stdout.split()
    assert not ignored, "README links to gitignored paths, which would 404 in a clone: " + ", ".join(ignored)


def test_public_container_dir_has_no_lab_tooling() -> None:
    """Host-specific job-distribution scripts are lab infrastructure and do not
    belong in a public repository."""
    names = {p.name for p in Path("container").iterdir() if p.is_file()}
    assert names == {
        "Dockerfile",
        "environment.container.yaml",
        "prepare-context.sh",
        "README.md",
    }
    for path in Path("container").glob("*"):
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="ignore")
            for leak in ["/Users/mtw", "venus", "mercury"]:
                assert leak not in text, f"{path.name} leaks '{leak}'"


def test_no_text_file_discloses_a_non_public_downstream_project() -> None:
    """The public repository must not describe what is built on top of it.

    Scope statements ("this project does not render HTML") are fine; naming a
    separate closed-source project, or characterising its business status, is
    not.
    """
    import re as _re

    disclosive = _re.compile(
        r"proprietary|closed[- ]source|commercial(ly)? licen|sibling project"
        r"|rnaconsnake[_-]portal",
        _re.IGNORECASE,
    )
    offenders: list[str] = []
    for path in _tracked_files():
        full = ROOT / path
        if not full.is_file():
            continue
        # This test names the forbidden words, so it cannot police itself.
        if full.resolve() == Path(__file__).resolve():
            continue
        try:
            text = full.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for match in disclosive.finditer(text):
            line = text[: match.start()].count("\n") + 1
            offenders.append(f"{path}:{line}: {match.group(0)!r}")
    assert not offenders, "public files disclose a non-public downstream project:\n" + "\n".join(offenders)


# --- refold ------------------------------------------------------------------
#
# A Python reimplementation of ViennaRNA's refold.pl. The constraint logic is
# pure text processing and is tested without the bindings; the folding tests
# skip where the bindings are absent.


def _refold_alignment(rows: dict[str, str]):
    from rnaconsnake.tools.alignment_io import Alignment

    return Alignment(order=list(rows), seqs=dict(rows))


def test_refold_pair_table_indexes_partners_and_rejects_unbalanced() -> None:
    from rnaconsnake.tools.refold import RefoldError, pair_table

    assert pair_table("((..))") == [5, 4, -1, -1, 1, 0]
    with pytest.raises(RefoldError):
        pair_table("((.)")
    with pytest.raises(RefoldError):
        pair_table(".)")


def test_refold_opens_pairs_to_a_gap_and_drops_the_column() -> None:
    """A pair whose partner is a gap in this sequence cannot form, and the gap
    column itself leaves the sequence, so both have to go."""
    from rnaconsnake.tools.refold import constrain

    sequence, constraint = constrain("GGC-AGCC", "(((())))", turn=0)
    assert sequence == "GGCAGCC"
    assert constraint == "(((.)))"


def test_refold_opens_pairs_the_sequence_cannot_form() -> None:
    from rnaconsnake.tools.refold import constrain

    # Columns 1 and 6 hold A and A, which is not a base pair, so that pair
    # opens while the enclosing G-C survives.
    sequence, constraint = constrain("GAGAAAAC", "((....))", turn=0)
    assert sequence == "GAGAAAAC"
    assert constraint == "(......)"

    # G-A on the outer pair too: nothing is left to constrain.
    _, constraint = constrain("GAGAAAAA", "((....))", turn=0)
    assert constraint == "........"


def test_refold_opens_pairs_closing_a_loop_shorter_than_turn() -> None:
    from rnaconsnake.tools.refold import constrain

    # The inner pair encloses two positions; the default turn of 3 opens it.
    _, tight = constrain("GGCCAAGGCC", "((((..))))", turn=3)
    assert tight == "(((....)))"
    _, loose = constrain("GGCCAAGGCC", "((((..))))", turn=2)
    assert loose == "((((..))))"


def test_refold_reads_the_consensus_from_a_dot_plot_above_the_threshold() -> None:
    from rnaconsnake.tools.refold import consensus_from_dotplot

    # The coloured layout RNAalifold -p writes, and the plain one.
    text = "\n".join(
        [
            "%!PS-Adobe-3.0 EPSF-3.0",
            "0.00 1.00 hsb 1 10 0.99 lbox",
            "0.00 1.00 hsb 2 9 0.50 lbox",
            "3 8 0.95 lbox",
            "1 10 0.99 ubox",
        ]
    )
    assert consensus_from_dotplot(text, 10, threshold=0.9) == "(.(....).)"
    assert consensus_from_dotplot(text, 10, threshold=0.4) == "(((....)))"


def test_refold_reads_the_consensus_from_the_rnaalifold_output() -> None:
    from rnaconsnake.tools.refold import RefoldError, consensus_from_alifold

    text = "\n".join(
        [
            ">cand_0001",
            "GGGCUAGCCC",
            "(((....))) (-5.20 = -4.10 + -1.10)",
            "(((....))) [-5.30]",
        ]
    )
    assert consensus_from_alifold(text) == "(((....)))"
    with pytest.raises(RefoldError):
        consensus_from_alifold("no structure here\n")


def test_refold_constraint_stream_matches_the_refold_pl_format() -> None:
    from rnaconsnake.tools.refold import format_constraints, refold_alignment

    alignment = _refold_alignment({"a": "GGCAAGCC", "b": "GGC-AGCC"})
    records = refold_alignment(alignment, "(((())))", turn=0)
    # a: the innermost A-A pair opens. b: the gap column goes and its partner
    # opens with it.
    assert format_constraints(records) == "> a\nGGCAAGCC\n(((..)))\n> b\nGGCAGCC\n(((.)))\n"


def test_refold_opens_a_pair_whose_partner_sits_in_column_zero() -> None:
    """refold.pl guards this branch with `$pt[$p] > 0`, leaves the `(` in
    column 0 unmatched, and then dies in its own pair-table check. Opening it
    is the only sensible reading, and it cannot change any output refold.pl
    produced successfully."""
    from rnaconsnake.tools.refold import constrain

    sequence, constraint = constrain("GAAA-", "(...)", turn=0)
    assert sequence == "GAAA"
    assert constraint == "...."


def test_refold_rejects_a_consensus_of_the_wrong_length() -> None:
    from rnaconsnake.tools.refold import RefoldError, constrain

    with pytest.raises(RefoldError):
        constrain("GGCC", "((((((", turn=0)


def test_refold_folds_under_the_constraint_and_formats_like_rnafold() -> None:
    """The folded output has to be what `refold.pl | RNAfold --noPS -C` writes,
    because that is the file the workflow records."""
    import re as _re

    pytest.importorskip("RNA")
    from rnaconsnake.tools.refold import fold_constrained, format_folded, refold_alignment

    alignment = _refold_alignment({"seq": "GGGGAAAACCCC"})
    records = refold_alignment(alignment, "((((....))))", turn=3)
    structure, energy = fold_constrained(records[0].sequence, records[0].constraint)
    assert structure == "((((....))))"
    assert energy < 0

    text = format_folded(records)
    lines = text.splitlines()
    assert lines[0] == "> seq"
    assert lines[1] == "GGGGAAAACCCC"
    # `structure ( -1.20)`: six columns of energy, as RNAfold prints it.
    assert _re.fullmatch(r"[().]{12} \( *-?\d+\.\d\d\)", lines[2]), lines[2]


def test_refold_matches_refold_pl_where_it_is_installed(tmp_path: Path) -> None:
    """Equivalence against the Perl original, when it is available to compare."""
    if shutil.which("refold.pl") is None or shutil.which("RNAfold") is None:
        pytest.skip("refold.pl and RNAfold are needed to compare against the original")
    pytest.importorskip("RNA")

    aln = tmp_path / "cand.aln"
    aln.write_text(
        "CLUSTAL 2.1 multiple sequence alignment\n\n"
        "seqA GGGCUAGCUAGGCAUCGAUCGGCUAGCUAGCCGAUCGAUGCCUAGCUAGCCC\n"
        "seqB GGGCUAGCUAGGCAUCGAUC-GCUAGCUAGCCGAUCGAUGCCUAGCUAGCCC\n"
        "seqC GGGCUAGCAAGGCAUCGAUCGGCUAGCUAGCCGAUCGAUGCCUAGCUUGCCC\n",
        encoding="utf-8",
    )
    consensus = tmp_path / "cand.alifold"
    consensus.write_text(
        ">cand\n"
        "GGGCUAGCUAGGCAUCGAUCGGCUAGCUAGCCGAUCGAUGCCUAGCUAGCCC\n"
        "((((((((((((((((((((((((....)))))))))))))))))))))))) (-58.97 = -49.63 + -9.35)\n",
        encoding="utf-8",
    )

    perl = subprocess.run(
        ["refold.pl", str(aln), str(consensus)],
        check=True,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
    )
    folded = subprocess.run(
        ["RNAfold", "--noPS", "-C"],
        input=perl.stdout,
        check=True,
        capture_output=True,
        text=True,
    )
    mine = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.refold",
            "--alignment",
            str(aln),
            "--consensus",
            str(consensus),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert mine.stdout == folded.stdout
