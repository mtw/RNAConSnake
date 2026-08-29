"""The export bundle: its tables, its locus and q-value
annotations, and the safety of the directory it is written to."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from rnaconsnake.workflow_helpers import CandidatePaths
from tests.helpers import (
    PYTHON,
    _write_nr_table,
    read_text,
    subprocess_env,
)


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
        paths.consensus_json: '{"alifold_consstruc": "<<>>"}\n',
        paths.refold_out: "> seqA\nACGU\n(..) ( -1.00)\n",
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
    # 1.2.0 exports the per-sequence refold and names the consensus artifact for
    # what it holds; 1.1.0 added the de-replication and calibration columns.
    assert manifest["export_schema_version"] == "1.2.0"
    assert manifest["dereplication"]["applied"] is False
    assert manifest["calibration"]["applied"] is False
    # Reported from the file the run was given, not assumed to be Stockholm.
    assert manifest["input_alignment_type"] == "stockholm"
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
    # The consensus artifact is named for what it holds, and the refold proper
    # -- computed on every candidate and previously never exported -- is there.
    assert "consensus_structure" in artifacts_text
    assert "Consensus secondary structure" in artifacts_text
    assert "refold_summary" not in artifacts_text
    assert f"{candidate_id}.refold.out" in artifacts_text
    assert "Per-sequence constrained refold" in artifacts_text
    assert read_text(
        bundle_dir / "files" / "np.mlocarna2020" / candidate_id / f"{candidate_id}.refold.out"
    ).startswith("> seqA")

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


def test_export_annotations_resolve_the_real_arm(tmp_path: Path) -> None:
    """With the null arm on, pipeline outputs live under arms/real/."""
    from rnaconsnake.export_annotations import analysis_root

    plain = tmp_path / "plain"
    (plain / "generated_files" / "summary" / "len_100").mkdir(parents=True)
    assert analysis_root(plain) == plain

    armed = tmp_path / "armed"
    (armed / "arms" / "real" / "generated_files" / "summary" / "len_100").mkdir(parents=True)
    assert analysis_root(armed) == armed / "arms" / "real"


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
    assert q is not None
    assert q.q_cascade == "0.1"
    assert q.cascade_pass == "yes"

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


def test_prepare_output_dir_refuses_a_directory_that_is_not_a_bundle(tmp_path: Path) -> None:
    """`RNAcs --export-bundle` used to force overwrite=True, so an existing
    directory was rmtree'd. Overwriting is now opt-in, and only ever replaces
    something that looks like a previous bundle."""
    from rnaconsnake.export_bundle import prepare_output_dir

    precious = tmp_path / "precious"
    precious.mkdir()
    (precious / "thesis.tex").write_text("do not delete\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="not part of an export bundle"):
        prepare_output_dir(precious, overwrite=True)
    assert (precious / "thesis.tex").is_file()

    with pytest.raises(FileExistsError, match="already exists and is not empty"):
        prepare_output_dir(precious, overwrite=False)
    assert (precious / "thesis.tex").is_file()


def test_prepare_output_dir_never_deletes_a_repository(tmp_path: Path) -> None:
    """The dangerous case: a stray dot entry must not be waved through as
    incidental. `.DS_Store` is; `.git` is emphatically not."""
    from rnaconsnake.export_bundle import prepare_output_dir

    repo = tmp_path / "repo"
    (repo / ".git").mkdir(parents=True)
    (repo / "manifest.json").write_text("{}\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match=r"\.git"):
        prepare_output_dir(repo, overwrite=True)
    assert (repo / ".git").is_dir()


def test_prepare_output_dir_ignores_a_finder_dotfile(tmp_path: Path) -> None:
    from rnaconsnake.export_bundle import prepare_output_dir

    bundle = tmp_path / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text("{}\n", encoding="utf-8")
    (bundle / ".DS_Store").write_bytes(b"\x00")

    prepare_output_dir(bundle, overwrite=True)

    assert list(bundle.iterdir()) == []


def test_prepare_output_dir_replaces_a_previous_bundle(tmp_path: Path) -> None:
    from rnaconsnake.export_bundle import prepare_output_dir

    bundle = tmp_path / "bundle"
    (bundle / "files").mkdir(parents=True)
    (bundle / "manifest.json").write_text("{}\n", encoding="utf-8")

    prepare_output_dir(bundle, overwrite=True)

    assert bundle.is_dir()
    assert list(bundle.iterdir()) == []


def test_prepare_output_dir_accepts_a_new_or_empty_directory(tmp_path: Path) -> None:
    from rnaconsnake.export_bundle import prepare_output_dir

    fresh = tmp_path / "fresh"
    prepare_output_dir(fresh, overwrite=False)
    assert fresh.is_dir()

    empty = tmp_path / "empty"
    empty.mkdir()
    prepare_output_dir(empty, overwrite=False)
    assert empty.is_dir()


def test_export_bundle_reports_the_input_alignment_format() -> None:
    """The manifest used to hard-code "stockholm", so a bundle from a Clustal
    run misstated its own provenance in the public contract."""
    from rnaconsnake.export_bundle import input_alignment_type

    assert input_alignment_type("/runs/aln.stk") == "stockholm"
    assert input_alignment_type("/runs/aln.STK") == "stockholm"
    assert input_alignment_type("/runs/aln.aln") == "clustal"
    assert input_alignment_type("/runs/aln.fasta") == "unknown"
    assert input_alignment_type(None) == "unknown"
