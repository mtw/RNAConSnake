"""Metric extraction from external-tool output, and the summary
log / CSV / Markdown reports built from it."""

from __future__ import annotations

import json
import subprocess
import textwrap
from pathlib import Path

from tests.helpers import (
    ALIFOLDZ_EMPTY_REPORT,
    ALIFOLDZ_GOOD_REPORT,
    PYTHON,
    _render_markdown,
    _run_extract_alifoldz,
    read_text,
    subprocess_env,
)


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


def test_extract_rnaz_reads_the_consensus_mfe_not_the_mean_single_sequence_mfe(tmp_path: Path) -> None:
    """RNAz prints "Mean single sequence MFE" *before* "Consensus MFE", and
    `re.search` returns the leftmost match. Matching both in one alternation
    recorded the mean single-sequence MFE as the consensus MFE on every real
    RNAz output -- a different quantity, in the summary tables and in the
    exported candidates.csv. The field ordering below is RNAz 2.1.1's."""
    rnaz = tmp_path / "one.rnaz.txt"
    out = tmp_path / "one.rnaz.json"
    rnaz.write_text(
        "\n".join(
            [
                "############################  RNAz 2.1.1  ##############################",
                "",
                " Sequences: 5",
                " Columns: 62",
                " Mean pairwise identity:  71.96",
                " Mean single sequence MFE: -20.93",
                " Consensus MFE: -18.74",
                " Energy contribution: -15.86",
                " Mean z-score:  -0.81",
                " Structure conservation index:   0.90",
                " SVM RNA-class probability: 0.916156",
                " Prediction: RNA",
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
    assert payload["consensus_mfe"] == "-18.74"
    assert payload["consensus_mfe"] != "-20.93"
    assert payload["sci"] == "0.90"
    assert payload["rnazprob"] == "0.916156"


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

    assert "rec_zero" in log_text
    assert "rscape 0" in log_text
    assert "rec_na" in log_text
    assert "rscape NA" in log_text
    assert "rec_zero,2,12,1,1,0,0.95,0.58,-14.10,-3.21,<<>>" in csv_text
    assert "rec_na,2,12,1,1,NA,0.95,0.58,-14.10,-3.21,<<>>" in csv_text
    assert "| rec_zero | 2 | 12 | 1 | 1 | 0 | 0.95 | 0.58 | -14.10 | -3.21 | <<>> |" in md_text
    assert "| rec_na | 2 | 12 | 1 | 1 | NA | 0.95 | 0.58 | -14.10 | -3.21 | <<>> |" in md_text
    assert "<<>>" in log_text


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
    header = next(line for line in text.splitlines() if line.startswith("| locus_id"))
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
    blocks = _re.split(r"^(?=\s*(?:rule|checkpoint)\s+[A-Za-z_])", text, flags=_re.MULTILINE)
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
