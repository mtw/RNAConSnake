from __future__ import annotations

import argparse
import csv
import re
import shutil
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from rnaconsnake import __version__
from rnaconsnake.export_annotations import analysis_root, manifest_block
from rnaconsnake.export_annotations import collect as collect_annotations
from rnaconsnake.workflow_helpers import SUMMARY_FIELDS, read_json, write_json

FEATURE_COLUMNS = [
    "feature_id",
    "feature_label",
    "feature_type",
    "candidate_count",
    "passing_count",
    "description",
    "source_label",
    "coordinate_label",
    "top_score",
    "covariation_supported_count",
    "rscape_count",
    "locus_count",
    "representative_count",
    "summary_md",
    "bundle_path",
    "last_updated",
]

CANDIDATE_COLUMNS = [
    "candidate_id",
    "feature_id",
    "start",
    "end",
    "length",
    "n_sequences",
    "status",
    "score_primary",
    "score_label",
    "rank_within_feature",
    "mean_pairwise_identity",
    "gap_fraction",
    "sci",
    "consensus_mfe",
    "covariation_pairs",
    "covariation_significant",
    "rscape_available",
    "cm_available",
    "rnazprob",
    "alifoldzscore",
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
    "summary_md",
    "notes",
]

ARTIFACT_COLUMNS = [
    "artifact_scope",
    "artifact_owner_id",
    "artifact_type",
    "artifact_label",
    "file_format",
    "path",
    "is_optional",
    "feature_id",
    "candidate_id",
    "preview_path",
    "group_name",
    "sort_order",
    "description",
]

CANDIDATE_COORDS = re.compile(r"_aln_(\d+)_(\d+)$")

CALIBRATION_ARTIFACTS = ["funnel.tsv", "qvalues.tsv", "score_dists.tsv", "summary.json"]


@dataclass(frozen=True)
class SummaryRecord:
    wlen: int
    values: dict[str, str]
    summary_path: Path

    @property
    def candidate_id(self) -> str:
        return self.values["wbn"]


def sanitize_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._-") or "dataset"


def write_csv(path: Path, columns: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def parse_float(value: str, default: float = float("-inf")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def discover_summary_records(run_dir: Path) -> list[SummaryRecord]:
    records: list[SummaryRecord] = []
    for path in sorted((run_dir / "generated_files" / "summary").glob("len_*/*.summary.json")):
        try:
            wlen = int(path.parent.name.split("_", 1)[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"Could not determine window length from {path}") from exc
        payload = read_json(path)
        normalized = {field: str(payload.get(field, "")) for field in SUMMARY_FIELDS}
        if not normalized["wbn"]:
            raise ValueError(f"Summary record missing wbn: {path}")
        records.append(SummaryRecord(wlen=wlen, values=normalized, summary_path=path))
    if not records:
        raise FileNotFoundError(
            f"No per-candidate summary JSON files found under {run_dir / 'generated_files' / 'summary'}"
        )
    return records


def sort_summary_records(records: list[SummaryRecord]) -> list[SummaryRecord]:
    """Same order as the run's own reports: strongest covariation first, then
    strongest AlifoldZ -- negated, because a more negative z-score is the more
    significant one. Missing values rank last."""
    return sorted(
        records,
        key=lambda record: (
            parse_float(record.values.get("maxcovarval", "")),
            -parse_float(record.values.get("alifoldzscore", ""), default=float("inf")),
        ),
        reverse=True,
    )


def copy_artifact(source: Path, bundle_root: Path, relative_path: str) -> str:
    destination = bundle_root / relative_path
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return relative_path.replace("\\", "/")


def detect_coords(candidate_id: str) -> tuple[int, int]:
    match = CANDIDATE_COORDS.search(candidate_id)
    if not match:
        raise ValueError(f"Could not extract start/end coordinates from candidate id: {candidate_id}")
    return int(match.group(1)), int(match.group(2))


def locus_columns(annotations, candidate_id: str) -> dict[str, object]:
    """De-replication and calibration columns for one candidate.

    Empty strings throughout when a run had neither, so a bundle from an
    uncalibrated run stays readable by the same consumer.
    """
    columns: dict[str, object] = {
        "locus_id": "",
        "locus_start": "",
        "locus_end": "",
        "locus_window_count": "",
        "is_representative": "",
        "redundant_to": "",
        "q_rnaz": "",
        "q_alifoldz": "",
        "q_cascade": "",
        "cascade_pass": "",
    }
    locus = annotations.locus_for(candidate_id)
    if locus is not None:
        representative = locus.is_representative(candidate_id)
        columns.update(
            {
                "locus_id": locus.locus_id,
                "locus_start": locus.locus_start,
                "locus_end": locus.locus_end,
                "locus_window_count": locus.window_count,
                "is_representative": bool_text(representative),
                # Empty for a representative; otherwise the candidate this one
                # is a fragment of, so nothing is orphaned.
                "redundant_to": "" if representative else locus.representative,
            }
        )
    qvalues = annotations.qvalues_for(candidate_id)
    if qvalues is not None:
        columns.update(
            {
                "q_rnaz": qvalues.q_rnaz,
                "q_alifoldz": qvalues.q_alifoldz,
                "q_cascade": qvalues.q_cascade,
                "cascade_pass": qvalues.cascade_pass,
            }
        )
    return columns


def bool_text(value: bool) -> str:
    return "true" if value else "false"


def write_methods_markdown(path: Path, feature_label: str, maxbpspans: list[int]) -> None:
    path.write_text(
        "\n".join(
            [
                "# RNAConSnake Export Bundle",
                "",
                f"- Feature label: {feature_label}",
                f"- Window lengths: {', '.join(str(wlen) for wlen in maxbpspans)}",
                "- Exported from RNAConSnake structured scientific outputs.",
                "- Visualization and browser-style rendering are intentionally outside "
                "the scope of this bundle.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_feature_markdown(
    path: Path, feature_row: dict[str, object], candidate_rows: list[dict[str, object]]
) -> None:
    path.write_text(
        "\n".join(
            [
                f"# Feature Summary: {feature_row['feature_label']}",
                "",
                f"- Feature ID: `{feature_row['feature_id']}`",
                f"- Candidate count: {feature_row['candidate_count']}",
                f"- Passing count: {feature_row['passing_count']}",
                f"- Covariation-supported count: {feature_row['covariation_supported_count']}",
                f"- R-scape count: {feature_row['rscape_count']}",
                "",
                "## Candidates",
                "",
            ]
            + [
                f"- `{row['candidate_id']}`: score `{row['score_primary']}`, status `{row['status']}`"
                for row in candidate_rows
            ]
            + [""]
        ),
        encoding="utf-8",
    )


def write_candidate_markdown(path: Path, record: SummaryRecord, candidate_row: dict[str, object]) -> None:
    values = record.values
    path.write_text(
        "\n".join(
            [
                f"# Candidate Summary: {candidate_row['candidate_id']}",
                "",
                f"- Feature ID: `{candidate_row['feature_id']}`",
                f"- Coordinates: {candidate_row['start']}..{candidate_row['end']}",
                f"- Length: {candidate_row['length']}",
                f"- Sequences: {candidate_row['n_sequences']}",
                f"- Status: `{candidate_row['status']}`",
                f"- Primary score (`{candidate_row['score_label']}`): {candidate_row['score_primary']}",
                f"- RNAz probability: {values['rnazprob']}",
                f"- alifoldz score: {values['alifoldzscore']}",
                f"- Max covariation value: {values['maxcovarval']}",
                f"- Max covariation count: {values['maxcovarcount']}",
                f"- R-scape covary count: {values['rscape_covary_count']}",
                "",
                "## Structures",
                "",
                f"- Consensus structure: `{values['alifold_consstruc']}`",
                "",
            ]
        ),
        encoding="utf-8",
    )


def add_artifact(
    artifacts: list[dict[str, object]],
    *,
    artifact_scope: str,
    artifact_owner_id: str,
    artifact_type: str,
    artifact_label: str,
    file_format: str,
    path: str,
    is_optional: bool,
    feature_id: str = "",
    candidate_id: str = "",
    preview_path: str = "",
    group_name: str = "",
    sort_order: int | str = "",
    description: str = "",
) -> None:
    artifacts.append(
        {
            "artifact_scope": artifact_scope,
            "artifact_owner_id": artifact_owner_id,
            "artifact_type": artifact_type,
            "artifact_label": artifact_label,
            "file_format": file_format,
            "path": path,
            "is_optional": bool_text(is_optional),
            "feature_id": feature_id,
            "candidate_id": candidate_id,
            "preview_path": preview_path,
            "group_name": group_name,
            "sort_order": sort_order,
            "description": description,
        }
    )


def build_export(args: argparse.Namespace) -> Path:
    run_dir = Path(args.run_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists():
        if not args.overwrite:
            raise FileExistsError(f"Output directory already exists: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    annotations = collect_annotations(run_dir)
    # With the null-model arm enabled the pipeline outputs live under
    # arms/real/; without it this is the run directory itself.
    summary_root = analysis_root(run_dir)
    records = sort_summary_records(discover_summary_records(summary_root))
    maxbpspans = sorted({record.wlen for record in records})

    feature_label_default = Path(args.input_alignment).stem if args.input_alignment else run_dir.name
    feature_id = sanitize_id(args.feature_id or feature_label_default)
    feature_label = args.feature_label or feature_label_default
    dataset_id = sanitize_id(args.dataset_id or run_dir.name)
    dataset_label = args.dataset_label or run_dir.name
    source_label = args.source_label or (Path(args.input_alignment).name if args.input_alignment else "")
    generated_at = datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    feature_summary_rel = f"feature_summaries/{feature_id}.md"
    methods_rel = "methods.md"
    candidate_rows: list[dict[str, object]] = []
    artifacts: list[dict[str, object]] = []

    for rank, record in enumerate(records, start=1):
        candidate_id = record.candidate_id
        start, end = detect_coords(candidate_id)
        candidate_bundle_dir = f"files/{feature_id}/{candidate_id}"
        candidate_summary_rel = f"candidate_summaries/{candidate_id}.md"

        maxcovar_value = parse_float(record.values["maxcovarval"], default=0.0)
        rscape_count = parse_int(record.values["rscape_covary_count"], default=0)

        candidate_row = {
            "candidate_id": candidate_id,
            "feature_id": feature_id,
            "start": start,
            "end": end,
            "length": end - start + 1,
            "n_sequences": parse_int(record.values["nrseq"], default=0),
            "status": "retained",
            "score_primary": f"{maxcovar_value:g}",
            "score_label": "maxcovarval",
            "rank_within_feature": rank,
            "mean_pairwise_identity": "",
            "gap_fraction": "",
            "sci": record.values.get("sci", ""),
            "consensus_mfe": record.values.get("consensus_mfe", ""),
            "covariation_pairs": str(rscape_count) if record.values["rscape_covary_count"] else "",
            "covariation_significant": bool_text(rscape_count > 0),
            "rscape_available": bool_text(False),
            "cm_available": bool_text(False),
            "rnazprob": record.values["rnazprob"],
            "alifoldzscore": record.values["alifoldzscore"],
            "summary_md": candidate_summary_rel,
            "notes": "",
        }
        candidate_row.update(locus_columns(annotations, candidate_id))

        source_files = {
            "stockholm_alignment": (
                summary_root / "Lalifold" / f"len_{record.wlen}" / "split" / f"{candidate_id}.stk",
                f"{candidate_bundle_dir}/{candidate_id}.split.stk",
                "Original split Stockholm alignment",
                "sto",
                False,
                "alignments",
                10,
            ),
            "filtered_alignment": (
                summary_root / "generated_files" / "stk" / f"len_{record.wlen}" / f"{candidate_id}.stk",
                f"{candidate_bundle_dir}/{candidate_id}.cleaned.stk",
                "Filtered cleaned Stockholm alignment",
                "sto",
                False,
                "alignments",
                20,
            ),
            "clustal_alignment": (
                summary_root / "generated_files" / "aln" / f"len_{record.wlen}" / f"{candidate_id}.aln",
                f"{candidate_bundle_dir}/{candidate_id}.clustal.aln",
                "Clustal reformatted alignment",
                "aln",
                False,
                "alignments",
                30,
            ),
            "metrics_table": (
                record.summary_path,
                f"{candidate_bundle_dir}/{candidate_id}.summary.json",
                "Structured candidate metrics",
                "json",
                False,
                "metrics",
                40,
            ),
            "refold_summary": (
                summary_root
                / "generated_files"
                / "refold"
                / f"len_{record.wlen}"
                / f"{candidate_id}.refold.json",
                f"{candidate_bundle_dir}/{candidate_id}.refold.json",
                "Refold summary",
                "json",
                False,
                "metrics",
                50,
            ),
            "alignment_plot_pdf": (
                summary_root
                / "generated_files"
                / "rnaalifold"
                / f"len_{record.wlen}"
                / candidate_id
                / f"{candidate_id}_aln.pdf",
                f"{candidate_bundle_dir}/{candidate_id}.alignment_plot.pdf",
                "Alignment plot PDF",
                "pdf",
                False,
                "figures",
                60,
            ),
            "consensus_plot_pdf": (
                summary_root
                / "generated_files"
                / "rnaalifold"
                / f"len_{record.wlen}"
                / candidate_id
                / f"{candidate_id}_ss.pdf",
                f"{candidate_bundle_dir}/{candidate_id}.consensus_plot.pdf",
                "Consensus secondary structure plot PDF",
                "pdf",
                False,
                "figures",
                70,
            ),
        }

        optional_files = {
            "rscape_pdf": (
                summary_root
                / "generated_files"
                / "rscape"
                / f"len_{record.wlen}"
                / f"{candidate_id}.sto.pdf",
                f"{candidate_bundle_dir}/{candidate_id}.rscape.pdf",
                "R-scape covariation PDF",
                "pdf",
                True,
                "covariation",
                80,
            ),
            "covariation_table": (
                summary_root / "generated_files" / "rscape" / f"len_{record.wlen}" / f"{candidate_id}.power",
                f"{candidate_bundle_dir}/{candidate_id}.rscape.power",
                "R-scape power output",
                "txt",
                True,
                "covariation",
                81,
            ),
        }

        for artifact_type, (
            source,
            relpath,
            label,
            file_format,
            is_optional,
            group_name,
            sort_order,
        ) in source_files.items():
            if not source.is_file():
                raise FileNotFoundError(f"Missing required export artifact for {candidate_id}: {source}")
            bundle_relpath = copy_artifact(source, output_dir, relpath)
            add_artifact(
                artifacts,
                artifact_scope="candidate",
                artifact_owner_id=candidate_id,
                artifact_type=artifact_type,
                artifact_label=label,
                file_format=file_format,
                path=bundle_relpath,
                is_optional=is_optional,
                feature_id=feature_id,
                candidate_id=candidate_id,
                group_name=group_name,
                sort_order=sort_order,
            )

        rscape_available = False
        for artifact_type, (
            source,
            relpath,
            label,
            file_format,
            is_optional,
            group_name,
            sort_order,
        ) in optional_files.items():
            if not source.is_file():
                continue
            bundle_relpath = copy_artifact(source, output_dir, relpath)
            add_artifact(
                artifacts,
                artifact_scope="candidate",
                artifact_owner_id=candidate_id,
                artifact_type=artifact_type,
                artifact_label=label,
                file_format=file_format,
                path=bundle_relpath,
                is_optional=is_optional,
                feature_id=feature_id,
                candidate_id=candidate_id,
                group_name=group_name,
                sort_order=sort_order,
            )
            rscape_available = True

        cm_status_path = (
            summary_root / "generated_files" / "cm" / f"len_{record.wlen}" / f"{candidate_id}.cm.status.json"
        )
        if cm_status_path.is_file():
            cm_status = read_json(cm_status_path)
            cm_relpath = str(cm_status.get("cm", "")).replace("\\", "/")
            if cm_relpath:
                cm_source = summary_root / cm_relpath
                if cm_source.is_file():
                    bundle_relpath = copy_artifact(
                        cm_source,
                        output_dir,
                        f"{candidate_bundle_dir}/{candidate_id}.cm",
                    )
                    add_artifact(
                        artifacts,
                        artifact_scope="candidate",
                        artifact_owner_id=candidate_id,
                        artifact_type="covariance_model",
                        artifact_label="Covariance model",
                        file_format="cm",
                        path=bundle_relpath,
                        is_optional=True,
                        feature_id=feature_id,
                        candidate_id=candidate_id,
                        group_name="models",
                        sort_order=90,
                    )
                    candidate_row["cm_available"] = bool_text(True)

        candidate_row["rscape_available"] = bool_text(rscape_available)

        candidate_md_path = output_dir / candidate_summary_rel
        candidate_md_path.parent.mkdir(parents=True, exist_ok=True)
        write_candidate_markdown(candidate_md_path, record, candidate_row)

        candidate_rows.append(candidate_row)

    feature_row = {
        "feature_id": feature_id,
        "feature_label": feature_label,
        "feature_type": args.feature_type,
        "candidate_count": len(candidate_rows),
        "passing_count": len(candidate_rows),
        "locus_count": annotations.locus_count,
        "representative_count": annotations.representative_count,
        "description": args.description or f"RNAConSnake export for {feature_label}",
        "source_label": source_label,
        "coordinate_label": "",
        "top_score": max(
            (parse_float(str(row["score_primary"]), default=0.0) for row in candidate_rows), default=0.0
        ),
        "covariation_supported_count": sum(
            1 for row in candidate_rows if parse_float(str(row["score_primary"]), default=0.0) > 0
        ),
        "rscape_count": sum(1 for row in candidate_rows if row["rscape_available"] == "true"),
        "summary_md": feature_summary_rel,
        "bundle_path": "",
        "last_updated": generated_at,
    }

    methods_path = output_dir / methods_rel
    write_methods_markdown(methods_path, feature_label=feature_label, maxbpspans=maxbpspans)
    feature_md_path = output_dir / feature_summary_rel
    feature_md_path.parent.mkdir(parents=True, exist_ok=True)
    write_feature_markdown(feature_md_path, feature_row, candidate_rows)

    for wlen in maxbpspans:
        for suffix, artifact_label, file_format, sort_order in [
            ("RNAConSnake.log", f"Summary log for len_{wlen}", "txt", 200),
            ("RNAConSnake.log.csv", f"Summary CSV for len_{wlen}", "csv", 201),
            ("RNAConSnake.md", f"Summary Markdown for len_{wlen}", "md", 202),
            ("RNAConSnake.nr.csv", f"Non-redundant loci for len_{wlen}", "csv", 203),
            ("RNAConSnake.nr.json", f"De-replication provenance for len_{wlen}", "json", 204),
        ]:
            source = summary_root / "generated_files" / "summary" / f"len_{wlen}" / suffix
            if source.is_file():
                bundle_relpath = copy_artifact(
                    source,
                    output_dir,
                    f"files/dataset/len_{wlen}/{suffix}",
                )
                add_artifact(
                    artifacts,
                    artifact_scope="dataset",
                    artifact_owner_id=dataset_id,
                    artifact_type="metrics_table",
                    artifact_label=artifact_label,
                    file_format=file_format,
                    path=bundle_relpath,
                    is_optional=False,
                    group_name="metrics",
                    sort_order=sort_order,
                )

    manifest = {
        "project": "RNAConSnake",
        "export_schema_version": "1.1.0",
        "pipeline_version": __version__,
        "dataset_id": dataset_id,
        "dataset_label": dataset_label,
        "generated_at": generated_at,
        "features_file": "features.csv",
        "candidates_file": "candidates.csv",
        "artifacts_file": "artifacts.csv",
        "description": args.description or f"RNAConSnake export for {feature_label}",
        "methods_file": methods_rel,
        "input_alignment_type": "stockholm",
        "feature_count": 1,
        "candidate_count": len(candidate_rows),
        "source_label": source_label,
        "notes": "",
        "features": {
            "methods_markdown": True,
            "feature_markdown": True,
            "candidate_markdown": True,
            "preview_images": False,
            "rscape_outputs": any(row["rscape_available"] == "true" for row in candidate_rows),
            "covariance_model_outputs": any(row["cm_available"] == "true" for row in candidate_rows),
            "dereplicated_loci": annotations.has_loci,
            "calibrated_qvalues": annotations.has_calibration,
        },
        **manifest_block(annotations),
    }

    write_json(output_dir / "manifest.json", manifest)
    if annotations.calibration_dir is not None:
        for index, name in enumerate(CALIBRATION_ARTIFACTS):
            source = annotations.calibration_dir / name
            if not source.is_file():
                continue
            bundle_relpath = copy_artifact(source, output_dir, f"files/dataset/calibration/{name}")
            add_artifact(
                artifacts,
                artifact_scope="dataset",
                artifact_owner_id=dataset_id,
                artifact_type="calibration_table",
                artifact_label=f"Null-model calibration: {name}",
                file_format=name.rsplit(".", 1)[-1],
                path=bundle_relpath,
                is_optional=True,
                group_name="calibration",
                sort_order=300 + index,
                description=(
                    "Empirical FDR / q-values from the null-model arm. "
                    "See the calibration block in manifest.json for whether the "
                    "reported FDR is conditional on stage-one survival."
                ),
            )

    if annotations.versions_file is not None:
        bundle_relpath = copy_artifact(annotations.versions_file, output_dir, "files/dataset/versions.yaml")
        add_artifact(
            artifacts,
            artifact_scope="dataset",
            artifact_owner_id=dataset_id,
            artifact_type="provenance",
            artifact_label="Toolchain versions",
            file_format="yaml",
            path=bundle_relpath,
            is_optional=True,
            group_name="provenance",
            sort_order=400,
            description="Exact external tool versions used for this run.",
        )

    write_csv(output_dir / "features.csv", FEATURE_COLUMNS, [feature_row])
    write_csv(output_dir / "candidates.csv", CANDIDATE_COLUMNS, candidate_rows)
    write_csv(output_dir / "artifacts.csv", ARTIFACT_COLUMNS, artifacts)

    validate_bundle(output_dir)
    return output_dir


def validate_bundle(bundle_root: Path) -> None:
    manifest = read_json(bundle_root / "manifest.json")
    features_path = bundle_root / manifest["features_file"]
    candidates_path = bundle_root / manifest["candidates_file"]
    artifacts_path = bundle_root / manifest["artifacts_file"]
    for path in [features_path, candidates_path, artifacts_path]:
        if not path.is_file():
            raise FileNotFoundError(f"Missing required bundle file: {path}")

    with features_path.open(encoding="utf-8", newline="") as handle:
        feature_rows = list(csv.DictReader(handle))
    with candidates_path.open(encoding="utf-8", newline="") as handle:
        candidate_rows = list(csv.DictReader(handle))
    with artifacts_path.open(encoding="utf-8", newline="") as handle:
        artifact_rows = list(csv.DictReader(handle))

    feature_ids = {row["feature_id"] for row in feature_rows}
    candidate_ids = {row["candidate_id"] for row in candidate_rows}
    if len(feature_ids) != len(feature_rows):
        raise ValueError("Feature IDs are not unique in features.csv")
    if len(candidate_ids) != len(candidate_rows):
        raise ValueError("Candidate IDs are not unique in candidates.csv")

    for row in candidate_rows:
        if row["feature_id"] not in feature_ids:
            raise ValueError(f"Candidate references unknown feature_id: {row['feature_id']}")

    for row in artifact_rows:
        artifact_path = bundle_root / row["path"]
        if not artifact_path.is_file():
            raise FileNotFoundError(f"Artifact path does not exist: {artifact_path}")
        preview_path = row.get("preview_path", "")
        if preview_path:
            resolved = bundle_root / preview_path
            if not resolved.is_file():
                raise FileNotFoundError(f"Preview path does not exist: {resolved}")
        scope = row["artifact_scope"]
        if scope == "feature" and row.get("feature_id", "") not in feature_ids:
            raise ValueError(
                f"Feature-scoped artifact references unknown feature_id: {row.get('feature_id', '')}"
            )
        if scope == "candidate" and row.get("candidate_id", "") not in candidate_ids:
            raise ValueError(
                f"Candidate-scoped artifact references unknown candidate_id: {row.get('candidate_id', '')}"
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Export a completed RNAConSnake run as a structured bundle.")
    parser.add_argument("run_dir", help="Completed RNAConSnake run directory")
    parser.add_argument("output_dir", help="Export bundle output directory")
    parser.add_argument("--dataset-id", help="Override machine-readable dataset identifier")
    parser.add_argument("--dataset-label", help="Override human-readable dataset label")
    parser.add_argument("--feature-id", help="Override feature identifier")
    parser.add_argument("--feature-label", help="Override feature display label")
    parser.add_argument(
        "--feature-type",
        default="other",
        help="Feature type label for features.csv, default: other",
    )
    parser.add_argument(
        "--input-alignment", help="Optional original input alignment path, used only for metadata defaults"
    )
    parser.add_argument("--source-label", help="Optional source/provenance label")
    parser.add_argument("--description", help="Optional dataset/feature description")
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite the output directory if it exists"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    build_export(args)
    print(f"RNAConSnake export bundle written to: {Path(args.output_dir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
