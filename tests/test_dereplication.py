"""Collapsing overlapping RNALalifold windows into loci, and
choosing the window that represents each one."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.helpers import (
    _candidate_rows,
    _window_rows,
    _write_summary_csv,
    read_text,
)


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
