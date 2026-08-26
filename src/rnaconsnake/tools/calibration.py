#!/usr/bin/env python3

"""Null-model calibration: empirical FDR / q-values and the filter funnel.

For each score ``s`` and threshold ``t``::

    FDR(t) = mean_over_replicates( N_null(s >= t) ) / N_real(s >= t)

and the reported q-value is the monotone envelope
``q(c) = min over t <= s(c) of FDR(t)``, clamped to ``[0, 1]``.

All counting happens on *merged loci* rather than raw RNALalifold windows (see
:mod:`rnaconsnake.tools.loci`), because overlapping windows are not
independent.  The collapse ratio is compared across arms and a warning is
recorded when the arms disagree by more than the configured tolerance: if
de-redundancy behaved differently on null alignments, the null count would
inflate and the FDR would be wrong in the flattering direction.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.tools.dereplicate import candidates_from_records, cluster_candidates
from rnaconsnake.tools.loci import collapse_ratio, parse_window_name
from rnaconsnake.workflow_helpers import REAL_ARM, arm_class, write_json


FUNNEL_STAGES = ["windows", "loci", "rnaz", "alifoldz", "rscape", "cascade"]

FUNNEL_COLUMNS = ["arm_class", "arm", "wlen", "stage", "n_in", "n_out"]

QVALUE_COLUMNS = [
    "locus_id",
    "wlen",
    "start",
    "end",
    "representative",
    "n_windows",
    "nrseq",
    "alilen",
    "rnazprob",
    "alifoldzscore",
    "rscape_covary_count",
    "q_rnaz",
    "q_alifoldz",
    "cascade_pass",
    "q_cascade",
]

SCORE_DIST_COLUMNS = ["arm", "wlen", "locus_id", "start", "end", "score", "value"]

NA = "NA"


@dataclass(frozen=True)
class Thresholds:
    rnaz_prob: float
    alifoldz: float
    rscape_min_pairs: int
    stage1_rnaz_prob: float
    locus_min_overlap: int
    collapse_ratio_tolerance: float
    dereplicate_method: str = "containment"
    pair_containment: float = 0.9
    max_container_width: int = 120
    container_min_coverage: float = 0.8
    representative_rule: str = "widest"

    def validate(self) -> None:
        if self.stage1_rnaz_prob > self.rnaz_prob:
            raise ValueError(
                "calibration.stage1_rnaz_prob must not exceed calibration.rnaz_prob_threshold, "
                f"got {self.stage1_rnaz_prob} > {self.rnaz_prob}: the reported RNAz filter would "
                "then select candidates whose AlifoldZ score was never computed."
            )


@dataclass(frozen=True)
class WindowRecord:
    arm: str
    wlen: int
    name: str
    start: int
    end: int
    rnazprob: float | None
    alifoldz: float | None
    rscape: int | None
    nrseq: str
    alilen: str


def _to_float(value: str | None) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return None if math.isnan(parsed) else parsed


def _to_int(value: str | None) -> int | None:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None


def read_summary_rows(path: str | Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8", newline="") as handle:
        return [row for row in csv.DictReader(handle) if (row.get("wbn") or "").strip()]


def window_record(row: dict[str, str], arm: str, wlen: int) -> WindowRecord:
    name = row["wbn"].strip()
    coords = parse_window_name(name)
    return WindowRecord(
        arm=arm,
        wlen=wlen,
        name=name,
        start=coords.start,
        end=coords.end,
        rnazprob=_to_float(row.get("rnazprob")),
        alifoldz=_to_float(row.get("alifoldzscore")),
        rscape=_to_int(row.get("rscape_covary_count")),
        nrseq=(row.get("nrseq") or "").strip(),
        alilen=(row.get("alilen") or "").strip(),
    )


@dataclass
class LocusRecord:
    arm: str
    wlen: int
    locus_id: str
    start: int
    end: int
    n_windows: int
    rep: WindowRecord

    def passes_rnaz(self, thresholds: Thresholds) -> bool:
        return self.rep.rnazprob is not None and self.rep.rnazprob >= thresholds.rnaz_prob

    def passes_alifoldz(self, thresholds: Thresholds) -> bool:
        return self.rep.alifoldz is not None and self.rep.alifoldz <= thresholds.alifoldz

    def passes_rscape(self, thresholds: Thresholds) -> bool:
        return self.rep.rscape is not None and self.rep.rscape >= thresholds.rscape_min_pairs

    def oriented(self, score: str) -> float | None:
        """Score oriented so that larger is always more significant."""
        if score == "rnaz_prob":
            return self.rep.rnazprob
        if score == "alifoldz":
            return None if self.rep.alifoldz is None else -self.rep.alifoldz
        raise ValueError(f"Unknown score {score!r}")


def build_loci(
    rows: list[dict[str, str]],
    arm: str,
    wlen: int,
    method: str,
    pair_threshold: float,
    min_overlap: int,
    max_container_width: int = 120,
    min_coverage: float = 0.8,
    representative_rule: str = "widest",
) -> list[LocusRecord]:
    """Cluster windows into loci using the *same* de-replication the ranked
    candidate table uses, so q-values describe the objects that get reported."""
    clusters = cluster_candidates(
        candidates_from_records(rows),
        method=method,
        pair_threshold=pair_threshold,
        min_overlap=min_overlap,
        label=f"len{wlen}",
        max_container_width=max_container_width,
        min_coverage=min_coverage,
        representative_rule=representative_rule,
    )
    return [
        LocusRecord(
            arm=arm,
            wlen=wlen,
            locus_id=cluster.locus_id,
            start=cluster.start,
            end=cluster.end,
            n_windows=len(cluster.members),
            rep=window_record(cluster.representative.record, arm, wlen),
        )
        for cluster in clusters
    ]


def rscape_was_evaluated(loci_by_arm: dict[str, list[LocusRecord]]) -> bool:
    """Did R-scape actually produce counts, or is the branch switched off?

    With ``do_rscape: false`` every candidate carries ``NA``, so keeping the
    R-scape filter in the cascade would zero the headline number for reasons
    that have nothing to do with the data.
    """
    return any(
        locus.rep.rscape is not None
        for loci in loci_by_arm.values()
        for locus in loci
    )


def funnel_counts(
    loci: list[LocusRecord],
    n_windows: int,
    thresholds: Thresholds,
    include_rscape: bool = True,
) -> dict[str, tuple[int, int]]:
    """Sequential cascade counts: ``stage -> (n_in, n_out)``."""
    n_loci = len(loci)
    after_rnaz = [locus for locus in loci if locus.passes_rnaz(thresholds)]
    after_alifoldz = [locus for locus in after_rnaz if locus.passes_alifoldz(thresholds)]
    after_rscape = (
        [locus for locus in after_alifoldz if locus.passes_rscape(thresholds)]
        if include_rscape
        else list(after_alifoldz)
    )
    return {
        "windows": (n_windows, n_windows),
        "loci": (n_windows, n_loci),
        "rnaz": (n_loci, len(after_rnaz)),
        "alifoldz": (len(after_rnaz), len(after_alifoldz)),
        "rscape": (len(after_alifoldz), len(after_rscape)),
        "cascade": (n_loci, len(after_rscape)),
    }


def cascade_survivors(
    loci: list[LocusRecord], thresholds: Thresholds, include_rscape: bool = True
) -> list[LocusRecord]:
    return [
        locus
        for locus in loci
        if locus.passes_rnaz(thresholds)
        and locus.passes_alifoldz(thresholds)
        and (not include_rscape or locus.passes_rscape(thresholds))
    ]


def empirical_fdr(
    real_scores: list[float],
    null_scores_per_replicate: list[list[float]],
) -> dict[float, float]:
    """``FDR(t)`` on the grid of observed real scores, clamped to ``[0, 1]``."""
    grid = sorted(set(real_scores))
    fdr: dict[float, float] = {}
    for threshold in grid:
        n_real = sum(1 for value in real_scores if value >= threshold)
        if n_real == 0:
            fdr[threshold] = 1.0
            continue
        null_counts = [
            sum(1 for value in replicate if value >= threshold)
            for replicate in null_scores_per_replicate
        ]
        mean_null = statistics.fmean(null_counts) if null_counts else 0.0
        fdr[threshold] = min(1.0, max(0.0, mean_null / n_real))
    return fdr


def qvalue_envelope(fdr: dict[float, float]) -> dict[float, float]:
    """Monotone-decreasing envelope: ``q(t) = min over t' <= t of FDR(t')``."""
    envelope: dict[float, float] = {}
    running = 1.0
    for threshold in sorted(fdr):
        running = min(running, fdr[threshold])
        envelope[threshold] = running
    return envelope


def _format_number(value: float | None, digits: int = 6) -> str:
    if value is None:
        return NA
    return f"{value:.{digits}g}"


def calibrate(
    arm_inputs: dict[str, dict[int, str]],
    thresholds: Thresholds,
    null_metadata: dict[str, object],
    output_dir: str | Path,
    two_stage: bool,
) -> dict[str, object]:
    thresholds.validate()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    arms = sorted(arm_inputs)
    if REAL_ARM not in arms:
        raise ValueError("Calibration requires a 'real' arm")
    null_arms = [arm for arm in arms if arm != REAL_ARM]
    if not null_arms:
        raise ValueError("Calibration requires at least one null replicate arm")

    loci_by_arm: dict[str, list[LocusRecord]] = {}
    windows_by_arm: dict[str, int] = {}
    loci_by_arm_wlen: dict[str, dict[int, list[LocusRecord]]] = {}
    windows_by_arm_wlen: dict[str, dict[int, int]] = {}

    for arm in arms:
        all_loci: list[LocusRecord] = []
        total_windows = 0
        loci_by_arm_wlen[arm] = {}
        windows_by_arm_wlen[arm] = {}
        for wlen in sorted(arm_inputs[arm]):
            rows = read_summary_rows(arm_inputs[arm][wlen])
            loci = build_loci(
                rows,
                arm,
                wlen,
                thresholds.dereplicate_method,
                thresholds.pair_containment,
                thresholds.locus_min_overlap,
                thresholds.max_container_width,
                thresholds.container_min_coverage,
                thresholds.representative_rule,
            )
            loci_by_arm_wlen[arm][wlen] = loci
            windows_by_arm_wlen[arm][wlen] = len(rows)
            all_loci.extend(loci)
            total_windows += len(rows)
        loci_by_arm[arm] = all_loci
        windows_by_arm[arm] = total_windows

    warnings: list[str] = []

    include_rscape = rscape_was_evaluated(loci_by_arm)
    if not include_rscape:
        warnings.append(
            "R-scape produced no covariation counts in any arm (do_rscape disabled, or no "
            "candidate was evaluated); the R-scape filter is reported but excluded from the "
            "cascade, which is therefore RNAz + AlifoldZ only"
        )

    funnel_by_arm: dict[str, dict[int, dict[str, tuple[int, int]]]] = {
        arm: {
            wlen: funnel_counts(
                loci_by_arm_wlen[arm][wlen],
                windows_by_arm_wlen[arm][wlen],
                thresholds,
                include_rscape=include_rscape,
            )
            for wlen in loci_by_arm_wlen[arm]
        }
        for arm in arms
    }

    # Non-independence guard: the window -> locus collapse must behave the same
    # way in both arms, otherwise the null count inflates and the FDR is wrong
    # in the flattering direction.
    real_ratio = collapse_ratio(windows_by_arm[REAL_ARM], len(loci_by_arm[REAL_ARM]))
    null_ratios = [
        collapse_ratio(windows_by_arm[arm], len(loci_by_arm[arm])) for arm in null_arms
    ]
    mean_null_ratio = statistics.fmean(null_ratios) if null_ratios else 0.0
    if real_ratio > 0:
        deviation = abs(mean_null_ratio / real_ratio - 1.0)
        if deviation > thresholds.collapse_ratio_tolerance:
            warnings.append(
                "window-to-locus collapse ratio differs between arms by "
                f"{deviation:.1%} (real={real_ratio:.3f}, null mean={mean_null_ratio:.3f}, "
                f"tolerance={thresholds.collapse_ratio_tolerance:.0%}); null survivor counts "
                "may be inflated and the reported FDR is not trustworthy"
            )
    else:
        warnings.append("real arm produced no loci; FDR is undefined")

    # --- per-score empirical FDR ------------------------------------------------
    envelopes: dict[str, dict[float, float]] = {}
    fdr_tables: dict[str, dict[float, float]] = {}
    for score in ["rnaz_prob", "alifoldz"]:
        real_values = [
            value
            for value in (locus.oriented(score) for locus in loci_by_arm[REAL_ARM])
            if value is not None
        ]
        null_values = [
            [
                value
                for value in (locus.oriented(score) for locus in loci_by_arm[arm])
                if value is not None
            ]
            for arm in null_arms
        ]
        fdr_tables[score] = empirical_fdr(real_values, null_values)
        envelopes[score] = qvalue_envelope(fdr_tables[score])

    # --- composite cascade ------------------------------------------------------
    real_cascade = cascade_survivors(loci_by_arm[REAL_ARM], thresholds, include_rscape)
    null_cascade_counts = [
        len(cascade_survivors(loci_by_arm[arm], thresholds, include_rscape))
        for arm in null_arms
    ]
    mean_null_cascade = statistics.fmean(null_cascade_counts) if null_cascade_counts else 0.0
    cascade_fdr: float | None
    if real_cascade:
        cascade_fdr = min(1.0, max(0.0, mean_null_cascade / len(real_cascade)))
    else:
        cascade_fdr = None
    cascade_ids = {locus.locus_id for locus in real_cascade}

    # --- outputs ----------------------------------------------------------------
    conditional = bool(two_stage)
    conditional_note = (
        "AlifoldZ was evaluated only on stage-one survivors "
        f"(RNAz class probability >= {thresholds.stage1_rnaz_prob}); the reported FDR is "
        "CONDITIONAL on passing stage one and must not be read as unconditional."
        if conditional
        else "AlifoldZ was evaluated on every candidate; the reported FDR is unconditional."
    )

    funnel_path = out / "funnel.tsv"
    with open(funnel_path, "w", encoding="utf-8", newline="") as handle:
        handle.write("# rnaconsnake calibration funnel\n")
        handle.write(f"# null_method\t{null_metadata.get('method', NA)}\n")
        handle.write(f"# replicates\t{len(null_arms)}\n")
        handle.write(f"# two_stage\t{str(conditional).lower()}\n")
        handle.write(f"# fdr_conditional_on_stage_one\t{str(conditional).lower()}\n")
        handle.write(f"# {conditional_note}\n")
        handle.write("# counts are on merged loci, not raw RNALalifold windows\n")
        handle.write(f"# rscape_in_cascade\t{str(include_rscape).lower()}\n")
        for message in warnings:
            handle.write(f"# WARNING\t{message}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(FUNNEL_COLUMNS)
        for arm in arms:
            for wlen in sorted(funnel_by_arm[arm]):
                for stage in FUNNEL_STAGES:
                    n_in, n_out = funnel_by_arm[arm][wlen][stage]
                    writer.writerow([arm_class(arm), arm, wlen, stage, n_in, n_out])
        all_wlens = sorted({wlen for arm in arms for wlen in funnel_by_arm[arm]})
        for wlen in all_wlens:
            for stage in FUNNEL_STAGES:
                ins = [funnel_by_arm[arm][wlen][stage][0] for arm in null_arms if wlen in funnel_by_arm[arm]]
                outs = [funnel_by_arm[arm][wlen][stage][1] for arm in null_arms if wlen in funnel_by_arm[arm]]
                writer.writerow(
                    [
                        "null_mean",
                        "null_mean",
                        wlen,
                        stage,
                        f"{statistics.fmean(ins):.4f}" if ins else NA,
                        f"{statistics.fmean(outs):.4f}" if outs else NA,
                    ]
                )

    qvalues_path = out / "qvalues.tsv"
    with open(qvalues_path, "w", encoding="utf-8", newline="") as handle:
        handle.write(f"# fdr_conditional_on_stage_one\t{str(conditional).lower()}\n")
        handle.write(f"# q_resolution\t{1.0 / len(null_arms):.6g}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(QVALUE_COLUMNS)
        for locus in sorted(loci_by_arm[REAL_ARM], key=lambda item: (item.wlen, item.start, item.end)):
            rnaz_score = locus.oriented("rnaz_prob")
            alifoldz_score = locus.oriented("alifoldz")
            q_rnaz = envelopes["rnaz_prob"].get(rnaz_score) if rnaz_score is not None else None
            q_alifoldz = (
                envelopes["alifoldz"].get(alifoldz_score) if alifoldz_score is not None else None
            )
            in_cascade = locus.locus_id in cascade_ids
            writer.writerow(
                [
                    locus.locus_id,
                    locus.wlen,
                    locus.start,
                    locus.end,
                    locus.rep.name,
                    locus.n_windows,
                    locus.rep.nrseq,
                    locus.rep.alilen,
                    _format_number(locus.rep.rnazprob),
                    _format_number(locus.rep.alifoldz),
                    NA if locus.rep.rscape is None else locus.rep.rscape,
                    _format_number(q_rnaz),
                    _format_number(q_alifoldz),
                    "yes" if in_cascade else "no",
                    _format_number(cascade_fdr) if in_cascade else NA,
                ]
            )

    dists_path = out / "score_dists.tsv"
    with open(dists_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SCORE_DIST_COLUMNS)
        for arm in null_arms:
            for locus in loci_by_arm[arm]:
                for score, raw in [
                    ("rnaz_prob", locus.rep.rnazprob),
                    ("alifoldz", locus.rep.alifoldz),
                    ("rscape_covary_count", locus.rep.rscape),
                ]:
                    if raw is None:
                        continue
                    writer.writerow(
                        [arm, locus.wlen, locus.locus_id, locus.start, locus.end, score, raw]
                    )

    summary = {
        "null": {
            "method": null_metadata.get("method"),
            "replicates": len(null_arms),
            "seed": null_metadata.get("seed"),
            "determinism": null_metadata.get("determinism"),
            "pool_sha256": null_metadata.get("pool_sha256"),
            "generator_warnings": null_metadata.get("warnings", []),
        },
        "two_stage": conditional,
        "fdr_conditional_on_stage_one": conditional,
        "fdr_conditionality_note": conditional_note,
        "thresholds": {
            "rnaz_prob": thresholds.rnaz_prob,
            "alifoldz": thresholds.alifoldz,
            "rscape_min_pairs": thresholds.rscape_min_pairs,
            "stage1_rnaz_prob": thresholds.stage1_rnaz_prob,
            "locus_min_overlap": thresholds.locus_min_overlap,
            "dereplicate_method": thresholds.dereplicate_method,
            "pair_containment": thresholds.pair_containment,
            # Clustering parameters change how many loci each arm reports, so
            # they change the q-values; a recorded calibration is only
            # reproducible with all of them.
            "max_container_width": thresholds.max_container_width,
            "container_min_coverage": thresholds.container_min_coverage,
            "representative_rule": thresholds.representative_rule,
            "collapse_ratio_tolerance": thresholds.collapse_ratio_tolerance,
        },
        "counting_unit": "merged_loci",
        "dereplication_method": thresholds.dereplicate_method,
        "representative_rule": thresholds.representative_rule,
        "rscape_evaluated": include_rscape,
        "cascade_filters": ["rnaz", "alifoldz"] + (["rscape"] if include_rscape else []),
        "q_resolution": 1.0 / len(null_arms),
        "collapse_ratio": {
            "real": real_ratio,
            "null_mean": mean_null_ratio,
            "tolerance": thresholds.collapse_ratio_tolerance,
        },
        "counts": {
            "real_windows": windows_by_arm[REAL_ARM],
            "real_loci": len(loci_by_arm[REAL_ARM]),
            "real_cascade_survivors": len(real_cascade),
            "null_cascade_survivors_mean": mean_null_cascade,
            "null_cascade_survivors_per_replicate": null_cascade_counts,
        },
        "cascade_fdr": cascade_fdr,
        "warnings": warnings,
        "outputs": {
            "funnel": str(funnel_path),
            "qvalues": str(qvalues_path),
            "score_dists": str(dists_path),
        },
    }
    write_json(out / "summary.json", summary)

    for message in warnings:
        print(f"WARNING: calibration: {message}", file=sys.stderr)
    return summary


def _parse_arm_input(token: str) -> tuple[str, int, str]:
    arm, wlen, path = token.split(":", 2)
    return arm, int(wlen), path


def main() -> int:
    parser = argparse.ArgumentParser(description="RNAConSnake null-model calibration.")
    parser.add_argument(
        "--arm-input",
        action="append",
        required=True,
        metavar="ARM:WLEN:PATH",
        help="Per-arm, per-window-length summary CSV. Repeat for every arm/window pair.",
    )
    parser.add_argument("--null-metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--rnaz-prob-threshold", type=float, required=True)
    parser.add_argument("--alifoldz-threshold", type=float, required=True)
    parser.add_argument("--rscape-min-pairs", type=int, default=1)
    parser.add_argument("--stage1-rnaz-prob", type=float, default=0.5)
    parser.add_argument("--locus-min-overlap", type=int, default=1)
    parser.add_argument("--dereplicate-method", default="containment")
    parser.add_argument("--max-container-width", type=int, default=120)
    parser.add_argument("--container-min-coverage", type=float, default=0.8)
    parser.add_argument("--representative", default="widest")
    parser.add_argument("--pair-containment", type=float, default=0.9)
    parser.add_argument("--collapse-ratio-tolerance", type=float, default=0.2)
    parser.add_argument("--two-stage", action="store_true")
    args = parser.parse_args()

    arm_inputs: dict[str, dict[int, str]] = {}
    for token in args.arm_input:
        arm, wlen, path = _parse_arm_input(token)
        arm_inputs.setdefault(arm, {})[wlen] = path

    null_metadata = json.loads(Path(args.null_metadata).read_text(encoding="utf-8"))

    calibrate(
        arm_inputs=arm_inputs,
        thresholds=Thresholds(
            rnaz_prob=args.rnaz_prob_threshold,
            alifoldz=args.alifoldz_threshold,
            rscape_min_pairs=args.rscape_min_pairs,
            stage1_rnaz_prob=args.stage1_rnaz_prob,
            locus_min_overlap=args.locus_min_overlap,
            collapse_ratio_tolerance=args.collapse_ratio_tolerance,
            dereplicate_method=args.dereplicate_method,
            pair_containment=args.pair_containment,
            max_container_width=args.max_container_width,
            container_min_coverage=args.container_min_coverage,
            representative_rule=args.representative,
        ),
        null_metadata=null_metadata,
        output_dir=args.output_dir,
        two_stage=args.two_stage,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
