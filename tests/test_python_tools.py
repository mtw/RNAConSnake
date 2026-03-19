from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from pathlib import Path


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


def test_cli_workflow_smoke_test_with_fake_rnalalifold(tmp_path: Path) -> None:
    data_dir = tmp_path / "data"
    data_dir.mkdir()
    (data_dir / "aln1.stk").write_text(
        "# STOCKHOLM 1.0\n#=GF ID dummy_input\nseq1 ACGU\n//\n",
        encoding="utf-8",
    )

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_rnalalifold(bin_dir)

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
            "--cores",
            "1",
            "--",
            "--directory",
            str(tmp_path),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=env,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "Lalifold" / "len_100" / "RC_100_0001.stk").is_file()
    assert (tmp_path / "Lalifold" / "len_200" / "RC_200_0001.stk").is_file()
    assert (tmp_path / "generated_files" / "remgap" / "len_100" / "manifest.txt").is_file()
    assert (tmp_path / "generated_files" / "strip" / "len_200" / "manifest.txt").is_file()
    assert read_text(tmp_path / "generated_files" / "strip" / "len_100" / "manifest.txt").strip() == (
        "RC_100_0001_aln_1_12_stripped.stk"
    )
