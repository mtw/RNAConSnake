"""The deterministic Stockholm helpers: parsing, splitting,
gap and redundancy filtering, case normalisation, and covariation counting."""

from __future__ import annotations

import subprocess
import textwrap
from pathlib import Path

import pytest

from tests.helpers import (
    FIXTURES,
    LOWERCASE_CLUSTAL,
    PYTHON,
    TOY_ALIGNMENT,
    read_text,
    subprocess_env,
)


def test_parse_stockholm_merges_gr_annotation_across_interleaved_blocks(tmp_path: Path) -> None:
    """``#=GR`` is per-column data, exactly like ``#=GC``. An interleaved record
    repeats it once per block; appending each block's line separately emitted
    several short ``#=GR`` rows against a full-length alignment, while the
    sequences and the ``#=GC`` line beside them had been merged correctly."""
    from rnaconsnake.tools.stockholm_utils import iter_stockholm_lines, parse_stockholm_records

    source = tmp_path / "interleaved.stk"
    source.write_text(
        "# STOCKHOLM 1.0\n"
        "seqA ACGU\n"
        "seqB ACGA\n"
        "#=GR seqA SS ....\n"
        "#=GR seqB SS ....\n"
        "#=GC SS_cons ....\n"
        "\n"
        "seqA ACGU\n"
        "seqB ACGA\n"
        "#=GR seqA SS (())\n"
        "#=GR seqB SS (())\n"
        "#=GC SS_cons (())\n"
        "//\n",
        encoding="utf-8",
    )

    (record,) = parse_stockholm_records(source)

    assert record.seqs == {"seqA": "ACGUACGU", "seqB": "ACGAACGA"}
    assert record.gc_lines == ["#=GC SS_cons ....(())"]
    assert record.gr_lines == [
        ("seqA", "#=GR seqA SS ....(())"),
        ("seqB", "#=GR seqB SS ....(())"),
    ]
    # Every annotation row spans the same columns as the alignment it annotates.
    lines = list(iter_stockholm_lines(record))
    widths = {len(line.split(None, 3)[3]) for line in lines if line.startswith("#=GR ")}
    assert widths == {8}


def test_parse_stockholm_keeps_distinct_gr_tags_apart(tmp_path: Path) -> None:
    from rnaconsnake.tools.stockholm_utils import parse_stockholm_records

    source = tmp_path / "tags.stk"
    source.write_text(
        "# STOCKHOLM 1.0\nseqA ACGU\n#=GR seqA SS ....\n#=GR seqA PP 9999\n//\n",
        encoding="utf-8",
    )

    (record,) = parse_stockholm_records(source)

    assert record.gr_lines == [("seqA", "#=GR seqA SS ...."), ("seqA", "#=GR seqA PP 9999")]


def test_parse_stockholm_still_rejects_a_malformed_gr_line(tmp_path: Path) -> None:
    from rnaconsnake.tools.stockholm_utils import parse_stockholm_records

    source = tmp_path / "bad.stk"
    source.write_text("# STOCKHOLM 1.0\nseqA ACGU\n#=GR seqA\n//\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Malformed #=GR line"):
        parse_stockholm_records(source)


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
    # Named for #=GF ID by default, even though this record also has an #=GF AC.
    assert read_text(tmp_path / "rec2.stk") == read_text(FIXTURES / "expected_split_acc2.stk")


def test_split_stockholm_accession_flag_selects_the_accession(tmp_path: Path) -> None:
    """``--accession`` used to decide nothing: both branches returned the
    accession, so an accession always won and the default never reached the
    ``#=GF ID``. RNALalifold writes no ``#=GF AC``, so the pipeline's own
    candidate names are the same either way."""
    for flags, expected, absent in [([], "rec2.stk", "acc2.stk"), (["--accession"], "acc2.stk", "rec2.stk")]:
        outdir = tmp_path / ("with_flag" if flags else "default")
        outdir.mkdir()
        subprocess.run(
            [
                PYTHON,
                "-m",
                "rnaconsnake.tools.split_stockholm",
                "-a",
                str(FIXTURES / "multi_input.stk"),
                *flags,
            ],
            cwd=outdir,
            check=True,
            env=subprocess_env(),
        )
        assert (outdir / expected).is_file()
        assert not (outdir / absent).exists()
        # The record with only an #=GF ID is named from it either way.
        assert (outdir / "rec1.stk").is_file()


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
    # The count is stdout, not the exit status. It used to be returned as the
    # status too, so a candidate with covariation exited non-zero and every
    # caller that checks -- run_checked, a `set -e` shell, CI -- read a
    # successful run as a failure.
    assert result.returncode == 0
    assert result.stdout == "3\tgapali_1"


def test_alifold_maxcovar_exits_zero_whatever_the_count() -> None:
    for stdin_text, expected in [
        ("  1 x 2 a b c d e\n  2 x 1 a b c d --:4\n", "3\tgapali_1"),
        ("  1 x 0 a b c d e f g h\n", "6\t1"),
        ("", "0\t0"),
    ]:
        result = subprocess.run(
            [PYTHON, "-m", "rnaconsnake.tools.alifold_maxcovar"],
            input=stdin_text,
            check=False,
            capture_output=True,
            text=True,
            env=subprocess_env(),
        )
        assert result.returncode == 0, f"exited {result.returncode} for {stdin_text!r}"
        assert result.stdout == expected


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
    assert "Seq_a" in result.stdout
    assert "SEQ_A" not in result.stdout
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
