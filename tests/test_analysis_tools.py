"""Analysis aids outside the workflow DAG: screenability
reporting and the sensitivity envelope."""

from __future__ import annotations

import subprocess
from pathlib import Path

from tests.helpers import (
    PYTHON,
    _envelope_alignment,
    read_text,
    subprocess_env,
)


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
    assert rows
    assert all(row["label"].startswith("n") for row in rows)

    subset = Subset(names=tuple(rows[0]["names"].split(";")), identity=rows[0]["identity"])
    path = tmp_path / "subset.stk"
    write_subset_alignment(alignment, subset, path, "demo")
    written = read_stockholm_alignment(path)
    assert written.length == alignment.length
    assert written.n_seq == subset.n_seq


def test_sensitivity_envelope_reports_having_no_subsets_instead_of_crashing(tmp_path: Path) -> None:
    """An alignment with fewer than two sequences has nothing to degrade. The
    final `min()` over no rows used to raise, after every output had been
    written."""
    alignment = tmp_path / "single.stk"
    alignment.write_text(
        "# STOCKHOLM 1.0\n#=GF ID lonely\na ACGUACGUACGU\n//\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "manifest.tsv"

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.sensitivity_envelope",
            "--alignment",
            str(alignment),
            "--output-dir",
            str(tmp_path / "subsets"),
            "--manifest",
            str(manifest),
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 1
    assert "no subsets to plan" in result.stderr
    assert "Traceback" not in result.stderr
    assert manifest.is_file()


def test_sensitivity_envelope_cli_writes_subsets_and_a_manifest(tmp_path: Path) -> None:
    alignment = tmp_path / "aln.stk"
    alignment.write_text(
        "# STOCKHOLM 1.0\n"
        "#=GF ID envelope\n"
        "a ACGUACGUACGUACGUACGU\n"
        "b ACGUACGUACGUACGUACGA\n"
        "c AGGUAGGUAGGUAGGUAGGU\n"
        "d UCAUUCAUUCAUUCAUUCAU\n"
        "//\n",
        encoding="utf-8",
    )
    outdir = tmp_path / "subsets"
    manifest = tmp_path / "manifest.tsv"

    result = subprocess.run(
        [
            PYTHON,
            "-m",
            "rnaconsnake.tools.sensitivity_envelope",
            "--alignment",
            str(alignment),
            "--output-dir",
            str(outdir),
            "--manifest",
            str(manifest),
            "--sizes",
            "2,3",
            "--identities",
            "0.4,0.8",
        ],
        check=False,
        capture_output=True,
        text=True,
        env=subprocess_env(),
    )

    assert result.returncode == 0, result.stderr
    assert "subsets spanning" in result.stdout
    text = read_text(manifest)
    assert text.startswith("# sensitivity envelope subsets of")
    rows = [line for line in text.splitlines() if not line.startswith("#")][1:]
    assert rows, "the manifest lists no subsets"
    # Every manifest row has an alignment written for it, with those sequences.
    from rnaconsnake.tools.alignment_io import read_stockholm_alignment

    for row in rows:
        label, n_seq, _identity, names = row.split("\t")
        subset = outdir / f"{label}.stk"
        assert subset.is_file(), label
        written = read_stockholm_alignment(subset)
        assert written.order == names.split(";")
        assert written.n_seq == int(n_seq)
        # Columns are preserved, so coordinates still mean the same thing.
        assert written.length == 20
