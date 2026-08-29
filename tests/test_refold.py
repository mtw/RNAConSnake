"""Constrained refolding, and folding an arbitrary alignment region."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from tests.helpers import (
    FOLD_REGION_ALIGNMENT,
    PYTHON,
    _fold_region_tools,
    _refold_alignment,
    read_text,
    subprocess_env,
    write_fake_postprocess_tools,
)


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


def test_slice_alignment_takes_inclusive_one_based_columns() -> None:
    from rnaconsnake.tools.alignment_io import Alignment
    from rnaconsnake.tools.fold_region import slice_alignment

    alignment = Alignment(order=["a", "b"], seqs={"a": "ACGUACGU", "b": "ACGU--GU"})
    region = slice_alignment(alignment, 2, 5)
    assert region.seqs == {"a": "CGUA", "b": "CGU-"}
    assert region.order == ["a", "b"]
    # The whole alignment is a legal span.
    assert slice_alignment(alignment, 1, 8).seqs == alignment.seqs


def test_slice_alignment_rejects_spans_outside_the_alignment() -> None:
    from rnaconsnake.tools.alignment_io import Alignment
    from rnaconsnake.tools.fold_region import slice_alignment

    alignment = Alignment(order=["a"], seqs={"a": "ACGUACGU"})
    for start, end in [(0, 4), (1, 9), (5, 4)]:
        with pytest.raises(ValueError, match="outside the alignment"):
            slice_alignment(alignment, start, end)


def test_write_stockholm_round_trips_through_the_reader(tmp_path: Path) -> None:
    from rnaconsnake.tools.alignment_io import Alignment, read_stockholm_alignment
    from rnaconsnake.tools.fold_region import write_stockholm

    alignment = Alignment(order=["a", "b"], seqs={"a": "ACGU", "b": "AC-U"})
    path = tmp_path / "region.stk"
    write_stockholm(alignment, path, "region_1_4")
    assert "#=GF ID region_1_4" in read_text(path)
    back = read_stockholm_alignment(path)
    assert back.order == alignment.order
    assert back.seqs == alignment.seqs


def test_fold_region_folds_plots_and_scores_a_span(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_postprocess_tools(bin_dir)
    alignment = tmp_path / "source.stk"
    alignment.write_text(FOLD_REGION_ALIGNMENT, encoding="utf-8")
    outdir = tmp_path / "out"

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.fold_region",
            "--alignment",
            str(alignment),
            "--start",
            "5",
            "--end",
            "24",
            "--output-dir",
            str(outdir),
            "--label",
            "elem",
            *_fold_region_tools(bin_dir),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stderr
    assert "elem: 5-24 (20 cols)" in result.stdout
    # The extracted span, the workflow's cleanup, the fold and both plots.
    assert (outdir / "elem.region.stk").is_file()
    assert (outdir / "elem.stk").is_file()
    assert (outdir / "elem.RNAalifold_results.stk").is_file()
    for produced in ["elem_aln.ps", "elem_ss.ps", "elem_aln.eps", "elem_aln.pdf", "elem_ss.pdf"]:
        assert (outdir / produced).is_file(), produced
    # The undeclared stray plot neither survives in outdir nor escapes into the
    # directory the tool was invoked from: alifoldz.pl drops one too.
    assert not (outdir / "alirna.ps").exists()
    assert not (Path.cwd() / "alirna.ps").exists()
    assert not (tmp_path / "alirna.ps").exists()

    payload = json.loads(read_text(outdir / "elem.scores.json"))
    assert payload["label"] == "elem"
    assert payload["start"] == 5
    assert payload["end"] == 24
    assert payload["columns"] == 20
    assert payload["n_seq_after_cleanup"] >= 1
    # Scores for this exact span, for comparison with what the screen reported.
    assert payload["rnazprob"] == "0.95"
    assert payload["alifoldzscore"] == "-3.21"


def test_fold_region_records_a_crashed_alifoldz_as_na(tmp_path: Path) -> None:
    """The 9999 sentinel is AlifoldZ's running-minimum initialiser, not a score."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_postprocess_tools(bin_dir)
    sentinel = bin_dir / "alifoldz.pl"
    sentinel.write_text(
        "#!/usr/bin/env python3\nimport sys\n_ = sys.stdin.read()\nprint('9999')\n",
        encoding="utf-8",
    )
    sentinel.chmod(0o755)

    alignment = tmp_path / "source.stk"
    alignment.write_text(FOLD_REGION_ALIGNMENT, encoding="utf-8")
    outdir = tmp_path / "out"
    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.fold_region",
            "--alignment",
            str(alignment),
            "--start",
            "1",
            "--end",
            "37",
            "--output-dir",
            str(outdir),
            "--label",
            "elem",
            *_fold_region_tools(bin_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    assert json.loads(read_text(outdir / "elem.scores.json"))["alifoldzscore"] == "NA"


def test_fold_region_skips_scoring_on_request(tmp_path: Path) -> None:
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_fake_postprocess_tools(bin_dir)
    alignment = tmp_path / "source.stk"
    alignment.write_text(FOLD_REGION_ALIGNMENT, encoding="utf-8")
    outdir = tmp_path / "out"

    subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.fold_region",
            "--alignment",
            str(alignment),
            "--start",
            "1",
            "--end",
            "37",
            "--output-dir",
            str(outdir),
            "--label",
            "elem",
            "--no-scores",
            *_fold_region_tools(bin_dir),
        ],
        check=True,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )
    payload = json.loads(read_text(outdir / "elem.scores.json"))
    assert "rnazprob" not in payload
    assert not (outdir / "elem.rnaz.txt").exists()
