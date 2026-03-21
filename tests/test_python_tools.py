from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path

from rnaconsnake import cli


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
            for i, arg in enumerate(args):
                if arg == "--id-prefix" and i + 1 < len(args):
                    prefix = args[i + 1]
                if arg == "-L" and i + 1 < len(args):
                    wlen = args[i + 1]

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
            sys.stdout.write(f"fake RNALalifold completed for window {wlen}\\n")
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
    assert result.stdout.strip() == "rnaconsnake 0.1.0"


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


def test_render_reports_writes_csv_markdown_and_html(tmp_path: Path) -> None:
    summary = tmp_path / "one.summary.json"
    rna_pdf_dir = tmp_path / "rnaalifold" / "len_100" / "rec1"
    rna_pdf_dir.mkdir(parents=True)
    (rna_pdf_dir / "rec1_aln.pdf").write_text("%PDF-FAKE\n", encoding="utf-8")
    (rna_pdf_dir / "rec1_ss.pdf").write_text("%PDF-FAKE\n", encoding="utf-8")
    rscape_pdf_dir = tmp_path / "rscape" / "len_100"
    rscape_pdf_dir.mkdir(parents=True)
    (rscape_pdf_dir / "rec1.sto.pdf").write_text("%PDF-FAKE\n", encoding="utf-8")
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
              "refold_firstseq": "ACGU",
              "refold_firststruc": "<<>>",
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
            "render-reports",
            "--label",
            "len_100",
            "--log",
            str(tmp_path / "pp_RNALalifold.log"),
            "--csv",
            str(tmp_path / "pp_RNALalifold.log.csv"),
            "--markdown",
            str(tmp_path / "pp_RNALalifold.md"),
            "--html",
            str(tmp_path / "summary" / "len_100" / "pp_RNALalifold.html"),
            str(summary),
        ],
        check=True,
        env=subprocess_env(),
    )

    assert "rec1" in read_text(tmp_path / "pp_RNALalifold.log")
    assert "wbn,nrseq,alilen" in read_text(tmp_path / "pp_RNALalifold.log.csv")
    assert "# RNAConSnake Summary: len_100" in read_text(tmp_path / "pp_RNALalifold.md")
    html_text = read_text(tmp_path / "summary" / "len_100" / "pp_RNALalifold.html")
    assert "<table>" in html_text
    assert "data-sort-type='number'" in html_text
    assert "Default sort: maxcovarval desc, then alifoldzscore desc." in html_text
    assert "filter-wbn" in html_text
    assert "filter-maxcovarval" in html_text
    assert "filter-rscape_covary_count" in html_text
    assert "plot-link" in html_text
    assert "aln pdf" in html_text
    assert "ss pdf" in html_text
    assert "data-preview-src='../../rnaalifold/len_100/rec1/rec1_aln.pdf'" in html_text
    assert "data-preview-src='../../rscape/len_100/rec1.sto.pdf'" in html_text
    assert "id='filter-alilen' type='number' step='any' value='20'" in html_text
    assert "id='filter-rnazprob' type='number' step='any' value='0.9'" in html_text
    assert "id='filter-alifoldzscore' type='number' step='any' value='-1'" in html_text
    assert "refold_firstseq" not in html_text
    assert "refold_firststruc" not in html_text
    assert "if (field === 'alifoldzscore')" in html_text
    assert "if (actual > threshold) return false;" in html_text


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


def test_cli_workflow_smoke_test_with_fake_rnalalifold(tmp_path: Path) -> None:
    input_alignment = tmp_path / "my_input.stk"
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
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "pp_RNALalifold.log").is_file()
    assert (tmp_path / "generated_files" / "summary" / "len_150" / "pp_RNALalifold.log.csv").is_file()
    assert (
        tmp_path
        / "generated_files"
        / "rnaalifold"
        / "len_150"
        / "RC_150_0001_aln_1_12"
        / "RC_150_0001_aln_1_12_aln.pdf"
    ).is_file()
    assert (
        tmp_path
        / "generated_files"
        / "rscape"
        / "len_150"
        / "RC_150_0001_aln_1_12.power"
    ).is_file()
    assert "RC_150_0001_aln_1_12" in read_text(
        tmp_path / "generated_files" / "summary" / "len_150" / "pp_RNALalifold.log"
    )
    assert "rscape_covary_count" in read_text(
        tmp_path / "generated_files" / "summary" / "len_150" / "pp_RNALalifold.log.csv"
    )
