#!/usr/bin/env python3

"""Locus and calibration annotations for the RNAConSnake export bundle.

The export bundle is RNAConSnake's public data contract: structured tables that
downstream consumers read. This module gathers the two pieces of information
that were previously visible only inside a run directory - which candidates are
de-replicated locus representatives, and what empirical q-values the null-model
calibration assigned - so that they travel with the bundle.

Everything here is *data*. Presentation of any kind belongs to whatever consumes
the bundle, never to RNAConSnake.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

CALIBRATION_DIRNAME = "calibration"
RESULTS_DIRNAME = "results"
REAL_ARM_DIRNAME = "arms/real"

CALIBRATION_FILES = ["funnel.tsv", "qvalues.tsv", "score_dists.tsv", "summary.json"]

NA_VALUES = {"", "NA", "na", "None", "nan"}


def analysis_root(run_dir: Path) -> Path:
    """Where the per-candidate outputs live.

    With the null-model arm enabled every pipeline output moves under
    ``arms/<arm>/``; the real arm is the one an export describes. Without it the
    run directory is the root, unchanged.
    """
    real_arm = run_dir / REAL_ARM_DIRNAME
    if (real_arm / "generated_files" / "summary").is_dir():
        return real_arm
    return run_dir


@dataclass(frozen=True)
class LocusAnnotation:
    """Which locus a window belongs to, after de-replication."""

    locus_id: str
    locus_start: int
    locus_end: int
    window_count: int
    representative: str
    members: tuple[str, ...] = ()

    def is_representative(self, candidate_id: str) -> bool:
        return candidate_id == self.representative


@dataclass(frozen=True)
class CalibrationAnnotation:
    """Empirical q-values assigned to one locus."""

    q_rnaz: str = ""
    q_alifoldz: str = ""
    q_cascade: str = ""
    cascade_pass: str = ""


@dataclass
class RunAnnotations:
    loci: dict[str, LocusAnnotation] = field(default_factory=dict)
    qvalues: dict[str, CalibrationAnnotation] = field(default_factory=dict)
    calibration_summary: dict = field(default_factory=dict)
    dereplication_method: str = ""
    calibration_dir: Path | None = None
    versions_file: Path | None = None
    nr_tables: dict[int, Path] = field(default_factory=dict)

    @property
    def has_loci(self) -> bool:
        return bool(self.loci)

    @property
    def has_calibration(self) -> bool:
        return bool(self.calibration_summary)

    @property
    def locus_count(self) -> int:
        return len({annotation.locus_id for annotation in self.loci.values()})

    @property
    def representative_count(self) -> int:
        return len({annotation.representative for annotation in self.loci.values()})

    def locus_for(self, candidate_id: str) -> LocusAnnotation | None:
        return self.loci.get(candidate_id)

    def qvalues_for(self, candidate_id: str) -> CalibrationAnnotation | None:
        annotation = self.loci.get(candidate_id)
        if annotation is None:
            return None
        return self.qvalues.get(annotation.locus_id)


def _clean(value: str | None) -> str:
    text = (value or "").strip()
    return "" if text in NA_VALUES else text


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def read_locus_annotations(
    root: Path,
) -> tuple[dict[str, LocusAnnotation], dict[int, Path], str]:
    """Map every window to its locus, from the per-window-length ``*.nr.csv``."""
    loci: dict[str, LocusAnnotation] = {}
    tables: dict[int, Path] = {}
    method = ""
    summary_dir = root / "generated_files" / "summary"
    for path in sorted(summary_dir.glob("len_*/RNAConSnake.nr.csv")):
        try:
            wlen = int(path.parent.name.split("_", 1)[1])
        except (IndexError, ValueError):
            continue
        tables[wlen] = path
        sidecar = path.with_suffix(".json")
        if sidecar.is_file():
            method = str(json.loads(sidecar.read_text(encoding="utf-8")).get("method", "")) or method
        with path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                representative = _clean(row.get("wbn"))
                if not representative:
                    continue
                members = tuple(member for member in (row.get("members") or "").split(";") if member)
                annotation = LocusAnnotation(
                    locus_id=_clean(row.get("locus_id")),
                    locus_start=int(row.get("locus_start") or 0),
                    locus_end=int(row.get("locus_end") or 0),
                    window_count=int(row.get("n_windows") or len(members) or 1),
                    representative=representative,
                    members=members or (representative,),
                )
                for member in annotation.members:
                    loci[member] = annotation
    return loci, tables, method


def read_calibration(run_dir: Path) -> tuple[dict[str, CalibrationAnnotation], dict, Path | None]:
    """Per-locus q-values plus the calibration summary, when the null arm ran."""
    calibration_dir = run_dir / RESULTS_DIRNAME / CALIBRATION_DIRNAME
    summary_path = calibration_dir / "summary.json"
    if not summary_path.is_file():
        return {}, {}, None

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    qvalues: dict[str, CalibrationAnnotation] = {}
    qvalues_path = calibration_dir / "qvalues.tsv"
    if qvalues_path.is_file():
        for row in _read_tsv(qvalues_path):
            locus_id = _clean(row.get("locus_id"))
            if not locus_id:
                continue
            qvalues[locus_id] = CalibrationAnnotation(
                q_rnaz=_clean(row.get("q_rnaz")),
                q_alifoldz=_clean(row.get("q_alifoldz")),
                q_cascade=_clean(row.get("q_cascade")),
                cascade_pass=_clean(row.get("cascade_pass")),
            )
    return qvalues, summary, calibration_dir


def collect(run_dir: Path) -> RunAnnotations:
    """Gather everything the bundle needs from a completed run directory."""
    run_dir = Path(run_dir)
    root = analysis_root(run_dir)
    loci, tables, method = read_locus_annotations(root)
    qvalues, summary, calibration_dir = read_calibration(run_dir)

    versions = run_dir / RESULTS_DIRNAME / "versions.yaml"
    return RunAnnotations(
        loci=loci,
        qvalues=qvalues,
        calibration_summary=summary,
        dereplication_method=method or str(summary.get("dereplication_method", "")),
        calibration_dir=calibration_dir,
        versions_file=versions if versions.is_file() else None,
        nr_tables=tables,
    )


def manifest_block(annotations: RunAnnotations) -> dict:
    """The ``dereplication`` / ``calibration`` manifest sections.

    Consumers must be able to tell, without reading any table, whether a bundle
    carries calibrated q-values and whether the reported FDR is conditional.
    """
    block: dict = {
        "dereplication": {
            "applied": annotations.has_loci,
            "method": annotations.dereplication_method,
            "locus_count": annotations.locus_count,
            "representative_count": annotations.representative_count,
        }
    }
    summary = annotations.calibration_summary
    if not summary:
        block["calibration"] = {"applied": False}
        return block

    null = summary.get("null", {}) or {}
    block["calibration"] = {
        "applied": True,
        "null_method": null.get("method", ""),
        "replicates": null.get("replicates", 0),
        "seed": null.get("seed", ""),
        "two_stage": bool(summary.get("two_stage", False)),
        "fdr_conditional_on_stage_one": bool(summary.get("fdr_conditional_on_stage_one", False)),
        "fdr_conditionality_note": summary.get("fdr_conditionality_note", ""),
        "cascade_fdr": summary.get("cascade_fdr"),
        "cascade_filters": summary.get("cascade_filters", []),
        "counting_unit": summary.get("counting_unit", ""),
        "q_resolution": summary.get("q_resolution"),
        "rscape_evaluated": summary.get("rscape_evaluated"),
        "thresholds": summary.get("thresholds", {}),
        "warnings": summary.get("warnings", []),
    }
    return block
