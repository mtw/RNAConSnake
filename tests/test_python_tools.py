from __future__ import annotations

import json
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from rnaconsnake import cli
from rnaconsnake.workflow_helpers import WorkflowSettings, CandidatePaths, initial_alignment_format_code


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
print("fake alifoldz header")
print("-3.21")
""",
        "RNAalifold": """\
#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
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
        "refold.pl": """\
#!/usr/bin/env python3
print(">fake")
print("ACGUACGU")
""",
        "RNAfold": """\
#!/usr/bin/env python3
import sys
_ = sys.stdin.read()
print(">fake")
print("ACGUACGU")
print("<<<<....>>>> (-1.0)")
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
    assert result.stdout.strip() == "rnaconsnake 0.2.0"


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
        "generated_files/rnaalifold/len_150/RC_150_0001_aln_1_12/"
        "RC_150_0001_aln_1_12.RNAalifold_results.stk"
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
    cm_status.write_text('{"built": true, "cm": "generated_files/cm/len_150/RC_150_0001_aln_1_12.cm"}\n', encoding="utf-8")
    (run_dir / "generated_files" / "cm" / "len_150" / f"{candidate_id}.cm").write_text("CM\n", encoding="utf-8")

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
    assert manifest["export_schema_version"] == "1.0.0"
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
    assert "artifact_scope,artifact_owner_id,artifact_type,artifact_label,file_format,path,is_optional" in artifacts_text
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
    assert (bundle_dir / "files" / "np.mlocarna2020" / candidate_id / f"{candidate_id}.alignment_plot.pdf").is_file()
    assert (bundle_dir / "files" / "np.mlocarna2020" / candidate_id / f"{candidate_id}.consensus_plot.pdf").is_file()


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
    assert csv_lines[1].startswith("rec_b,")
    assert csv_lines[2].startswith("rec_c,")
    assert csv_lines[3].startswith("rec_a,")


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

    assert read_text(tmp_path / "RC_100_0001_aln_2294_2375.stk") == read_text(real / "RC_100_0001_aln_2294_2375.stk")
    assert read_text(tmp_path / "RC_100_0001_aln_4042_4067.stk") == read_text(real / "RC_100_0001_aln_4042_4067.stk")


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


def test_failed_alifoldz_fallback_uses_numeric_zero(tmp_path: Path) -> None:
    out = tmp_path / "failed.alifoldz.json"
    out.write_text('{"alifoldzscore": "0.0"}\n', encoding="utf-8")

    assert '"alifoldzscore": "0.0"' in read_text(out)


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
    assert not (tmp_path / "generated_files" / "summary" / "len_150" / "RNAConSnake.html").exists()
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
    assert (
        tmp_path
        / "generated_files"
        / "rscape"
        / "len_150"
        / "RC_150_0001_aln_1_12.power"
    ).is_file()
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
