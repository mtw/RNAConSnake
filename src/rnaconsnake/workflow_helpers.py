from __future__ import annotations

import contextlib
import hashlib
import json
import shlex
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SUMMARY_FIELDS = [
    "wbn",
    "nrseq",
    "alilen",
    "maxcovarval",
    "maxcovarcount",
    "rscape_covary_count",
    "rscape_avg_confidence",
    "rscape_mutual_info",
    "rnazprob",
    "sci",
    "consensus_mfe",
    "alifoldzscore",
    "alifold_consstruc",
]


@dataclass(frozen=True)
class WorkflowSettings:
    input_alignment: str | None
    maxbpspan: list[int]
    lalifold_threads: int
    remove_gaponly_gapratio: float
    remove_gaponly_max_n: int
    do_cm: bool
    do_locarnate: bool
    do_png: bool
    do_rscape: bool
    rnaz_no_shuffle: bool
    cm_rnaz_threshold: float
    cm_alifoldz_threshold: float
    tools: dict[str, str]

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> WorkflowSettings:
        maxbpspan_vals = [int(value) for value in config.get("maxbpspan", [100, 200])]
        for val in maxbpspan_vals:
            if val <= 0 or val > 10000:
                raise ValueError(f"--maxbpspan {val} out of bounds: must be 1-10000 nt")
        return cls(
            input_alignment=config.get("input_alignment"),
            maxbpspan=maxbpspan_vals,
            lalifold_threads=int(config.get("lalifold_threads", 1)),
            remove_gaponly_gapratio=float(config.get("remove_gaponly_gapratio", 0.5)),
            remove_gaponly_max_n=int(config.get("remove_gaponly_max_n", 0)),
            do_cm=bool(config.get("do_cm", False)),
            do_locarnate=bool(config.get("do_locarnate", False)),
            do_png=bool(config.get("do_png", True)),
            do_rscape=bool(config.get("do_rscape", False)),
            rnaz_no_shuffle=bool(config.get("rnaz_no_shuffle", True)),
            cm_rnaz_threshold=float(config.get("cm_rnaz_prob_threshold", 0.9)),
            cm_alifoldz_threshold=float(config.get("cm_alifoldz_threshold", -2.0)),
            tools={str(key): str(value) for key, value in dict(config.get("tools", {})).items()},
        )

    def command_tokens(self, name: str, default: str) -> list[str]:
        return shlex.split(str(self.tools.get(name, default)))


NULL_METHODS = ("sissiz", "rnazRandomizeAln", "none")

REAL_ARM = "real"

ARM_WILDCARD_PATTERN = r"real|null_\d{3}"


@dataclass(frozen=True)
class NullSettings:
    """Configuration of the null-model calibration arm.

    The null arm is *off* by default: ``method: none`` (or ``replicates: 0``)
    reproduces the pre-calibration pipeline byte-for-byte, because the whole
    ``arms/{arm}/`` path prefix collapses to the empty string.
    """

    method: str
    replicates: int
    seed: int
    two_stage: bool
    pool_file: str | None = None

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> NullSettings:
        # An unquoted ``null:`` key in YAML parses as the null scalar, so the
        # section can arrive under ``None`` instead of ``"null"``. Accept both
        # rather than silently ignoring a user's configuration.
        section = config.get("null")
        if section is None:
            section = config.get(None)
        raw = dict(section or {})
        method = str(raw.get("method", "none"))
        if method not in NULL_METHODS:
            raise ValueError(f"Unknown null.method {method!r}. Expected one of: " + ", ".join(NULL_METHODS))
        replicates = int(raw.get("replicates", 10))
        if replicates > 1000:
            import warnings
            msg = (
                f"--null-replicates {replicates} is very large; "
                "this may take a very long time to run"
            )
            warnings.warn(msg, UserWarning, stacklevel=2)
        return cls(
            method=method,
            replicates=replicates,
            seed=int(raw.get("seed", 20261101)),
            two_stage=bool(raw.get("two_stage", True)),
            pool_file=_optional_path(raw.get("pool_file")),
        )

    @property
    def enabled(self) -> bool:
        return self.method != "none" and self.replicates > 0

    @property
    def effective_replicates(self) -> int:
        return self.replicates if self.enabled else 0

    def null_arms(self) -> list[str]:
        return [f"null_{index:03d}" for index in range(self.effective_replicates)]

    def arms(self) -> list[str]:
        if not self.enabled:
            return []
        return [REAL_ARM, *self.null_arms()]

    def arm_seed(self, arm: str) -> int | None:
        return arm_seed(arm, self.seed)


@dataclass(frozen=True)
class CalibrationSettings:
    """Thresholds and tolerances used by the calibration aggregation step."""

    rnaz_prob_threshold: float
    alifoldz_threshold: float
    rscape_min_pairs: int
    rscape_min_confidence: float
    rscape_min_mutual_info: float
    stage1_rnaz_prob: float
    locus_min_overlap: int
    collapse_ratio_tolerance: float
    dereplicate_method: str
    pair_containment: float
    max_container_width: int
    container_min_coverage: float
    representative_rule: str

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> CalibrationSettings:
        raw = dict(config.get("calibration", {}) or {})
        dereplicate = dict(config.get("dereplicate", {}) or {})
        return cls(
            rnaz_prob_threshold=float(
                raw.get("rnaz_prob_threshold", config.get("cm_rnaz_prob_threshold", 0.9))
            ),
            alifoldz_threshold=float(
                raw.get("alifoldz_threshold", config.get("cm_alifoldz_threshold", -2.0))
            ),
            rscape_min_pairs=int(raw.get("rscape_min_pairs", 1)),
            rscape_min_confidence=float(raw.get("rscape_min_confidence", 0.5)),
            rscape_min_mutual_info=float(raw.get("rscape_min_mutual_info", 0.1)),
            stage1_rnaz_prob=float(raw.get("stage1_rnaz_prob", 0.5)),
            locus_min_overlap=int(raw.get("locus_min_overlap", 1)),
            collapse_ratio_tolerance=float(raw.get("collapse_ratio_tolerance", 0.2)),
            dereplicate_method=str(dereplicate.get("method", "containment")),
            pair_containment=float(dereplicate.get("pair_containment", 0.9)),
            max_container_width=int(dereplicate.get("max_container_width", 120)),
            container_min_coverage=float(dereplicate.get("container_min_coverage", 0.8)),
            representative_rule=str(dereplicate.get("representative", "widest")),
        )

    def as_dict(self) -> dict[str, Any]:
        return {
            "rnaz_prob_threshold": self.rnaz_prob_threshold,
            "alifoldz_threshold": self.alifoldz_threshold,
            "rscape_min_pairs": self.rscape_min_pairs,
            "rscape_min_confidence": self.rscape_min_confidence,
            "rscape_min_mutual_info": self.rscape_min_mutual_info,
            "stage1_rnaz_prob": self.stage1_rnaz_prob,
            "locus_min_overlap": self.locus_min_overlap,
            "collapse_ratio_tolerance": self.collapse_ratio_tolerance,
            "dereplicate_method": self.dereplicate_method,
            "pair_containment": self.pair_containment,
            "max_container_width": self.max_container_width,
            "container_min_coverage": self.container_min_coverage,
            "representative_rule": self.representative_rule,
        }


def _optional_path(value: Any) -> str | None:
    """Normalise an optional path from config.

    Snakemake stringifies nested ``--config`` values, so a ``None`` arrives as
    the literal string ``"None"``, which is truthy. Treat the stringified
    empties as absent rather than as a filename.
    """
    if value is None:
        return None
    text = str(value).strip()
    return None if text in {"", "None", "null", "~"} else text


def arm_seed(arm: str, base_seed: int) -> int | None:
    """Deterministic per-arm seed derived from the arm name and the base seed.

    Returns ``None`` for the real arm, which is never simulated.
    """
    if arm == REAL_ARM:
        return None
    digest = hashlib.sha256(f"{int(base_seed)}:{arm}".encode()).hexdigest()
    return int(digest[:8], 16)


# Seeds Perl's RNG before handing control to a helper script, so that a rerun
# reproduces the same shuffles. Used for alifoldz.pl and rnazRandomizeAln.pl,
# neither of which exposes a seed option.
PERL_SEED_BOOTSTRAP = "srand(shift); my $script = shift; do $script; die $@ if $@;"


def perl_seed_env(env: dict[str, str] | None = None) -> dict[str, str]:
    """Environment that makes a seeded Perl run reproducible.

    Perl randomises hash iteration order per process, so ``srand`` alone is not
    enough to pin behaviour that walks a hash.
    """
    import os

    resolved = dict(env if env is not None else os.environ)
    resolved["PERL_HASH_SEED"] = "0"
    resolved["PERL_PERTURB_KEYS"] = "0"
    return resolved


def perl_seeded_command(
    tokens: list[str], seed: int, args: list[str] | None = None
) -> tuple[list[str], bool]:
    """Build a seeded invocation of a Perl helper script.

    Returns ``(command, seeded)``. Falls back to the command unchanged, and
    reports ``seeded=False``, when it is not a Perl script -- CI fakes and site
    wrappers keep working, they just are not reproducible.
    """
    import shutil as _shutil

    if not tokens:
        raise ValueError("Empty command; cannot seed")
    if Path(tokens[0]).name in {"perl", "perl5"} and len(tokens) > 1:
        script, extra = tokens[1], tokens[2:]
    else:
        script = _shutil.which(tokens[0]) or tokens[0]
        extra = tokens[1:]
    if not script.endswith(".pl"):
        return [*tokens, *(args or [])], False
    return ["perl", "-e", PERL_SEED_BOOTSTRAP, str(seed), script, *extra, *(args or [])], True


def derived_seed(base_seed: int, label: str) -> int:
    """Deterministic per-item seed, so parallel jobs differ but reproduce."""
    digest = hashlib.sha256(f"{int(base_seed)}:{label}".encode()).hexdigest()
    return int(digest[:8], 16)


def arm_prefix_for(arm: str | None) -> str:
    """Path prefix for an arm; the empty string when the null arm is disabled."""
    if not arm:
        return ""
    return f"arms/{arm}/"


def arm_class(arm: str) -> str:
    return REAL_ARM if arm == REAL_ARM else "null"


# The output roots the pipeline writes beneath. Finding any one of them under
# ``arms/real/`` settles which layout a finished run used.
ANALYSIS_ROOT_MARKERS = ("generated_files", "Lalifold")


def analysis_root(run_dir: str | Path) -> Path:
    """Where a completed run's per-candidate outputs live.

    With the null-model arm enabled every pipeline output moves under
    ``arms/<arm>/``, and the real arm is the one a report or an export
    describes. Without it the run directory is the root, unchanged. Anything
    reading a finished run has to resolve this first: looking straight at the
    run directory finds nothing at all in a calibrated run.
    """
    run_dir = Path(run_dir)
    real_arm = run_dir / "arms" / REAL_ARM
    if any((real_arm / marker).is_dir() for marker in ANALYSIS_ROOT_MARKERS):
        return real_arm
    return run_dir


@dataclass(frozen=True)
class CandidatePaths:
    wlen: str | int
    file: str
    arm_prefix: str = ""

    @property
    def len_dir(self) -> str:
        return f"len_{self.wlen}"

    @property
    def arm_alignment(self) -> str:
        return f"{self.arm_prefix}input_alignment.stk"

    @property
    def split(self) -> str:
        return f"{self.arm_prefix}Lalifold/{self.len_dir}/split/{self.file}.stk"

    @property
    def orig(self) -> str:
        return f"{self.arm_prefix}generated_files/orig/{self.len_dir}/{self.file}.orig.stk"

    @property
    def remgap(self) -> str:
        return f"{self.arm_prefix}generated_files/remgap/{self.len_dir}/{self.file}_remgap.stk"

    @property
    def strip(self) -> str:
        return f"{self.arm_prefix}generated_files/strip/{self.len_dir}/{self.file}_stripped.stk"

    @property
    def stk(self) -> str:
        return f"{self.arm_prefix}generated_files/stk/{self.len_dir}/{self.file}.stk"

    @property
    def aln(self) -> str:
        return f"{self.arm_prefix}generated_files/aln/{self.len_dir}/{self.file}.aln"

    @property
    def rnaz_txt(self) -> str:
        return f"{self.arm_prefix}generated_files/rnaz/{self.len_dir}/{self.file}.rnaz.txt"

    @property
    def rnaz_json(self) -> str:
        return f"{self.arm_prefix}generated_files/rnaz/{self.len_dir}/{self.file}.rnaz.json"

    @property
    def alifoldz_txt(self) -> str:
        return f"{self.arm_prefix}generated_files/alifoldz/{self.len_dir}/{self.file}.alifoldz.txt"

    @property
    def alifoldz_json(self) -> str:
        return f"{self.arm_prefix}generated_files/alifoldz/{self.len_dir}/{self.file}.alifoldz.json"

    @property
    def rnaalifold_dir(self) -> str:
        return f"{self.arm_prefix}generated_files/rnaalifold/{self.len_dir}/{self.file}"

    @property
    def rnaalifold_stdout(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}.alifold.out"

    @property
    def rnaalifold_stk(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}.RNAalifold_results.stk"

    @property
    def rnaalifold_aln(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}.RNAalifold_results.aln"

    @property
    def rnaalifold_aln_backup(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}.RNAalifold_results.aln~"

    @property
    def rnaalifold_cleaned_aln(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}.RNAalifold_results.cleaned.aln"

    @property
    def ali_out(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_ali.out"

    @property
    def dp_ps(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_dp.ps"

    @property
    def aln_ps(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_aln.ps"

    @property
    def aln_eps(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_aln.eps"

    @property
    def aln_pdf(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_aln.pdf"

    @property
    def ss_ps(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_ss.ps"

    @property
    def ss_eps(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_ss.eps"

    @property
    def ss_pdf(self) -> str:
        return f"{self.rnaalifold_dir}/{self.file}_ss.pdf"

    @property
    def refold_out(self) -> str:
        return f"{self.arm_prefix}generated_files/refold/{self.len_dir}/{self.file}_refold.out"

    @property
    def consensus_json(self) -> str:
        # The RNAalifold consensus structure, not the refold. It lived at
        # `refold/{file}.refold.json` and was named for the leg that runs beside
        # it rather than for what it holds; `_refold.out` is the refold.
        return f"{self.arm_prefix}generated_files/consensus/{self.len_dir}/{self.file}.consensus.json"

    @property
    def maxcovar_log(self) -> str:
        return f"{self.arm_prefix}generated_files/maxcovar/{self.len_dir}/{self.file}_alifoldmaxcovar.log"

    @property
    def maxcovar_json(self) -> str:
        return f"{self.arm_prefix}generated_files/maxcovar/{self.len_dir}/{self.file}.maxcovar.json"

    @property
    def rscape_power(self) -> str:
        return f"{self.arm_prefix}generated_files/rscape/{self.len_dir}/{self.file}.power"

    @property
    def rscape_json(self) -> str:
        return f"{self.arm_prefix}generated_files/rscape/{self.len_dir}/{self.file}.rscape.json"

    @property
    def rscape_sto_pdf(self) -> str:
        return f"{self.arm_prefix}generated_files/rscape/{self.len_dir}/{self.file}.sto.pdf"

    @property
    def summary_json(self) -> str:
        return f"{self.arm_prefix}generated_files/summary/{self.len_dir}/{self.file}.summary.json"

    @property
    def png_aln(self) -> str:
        return f"{self.arm_prefix}generated_files/png/{self.len_dir}/{self.file}_aln.png"

    @property
    def png_ss(self) -> str:
        return f"{self.arm_prefix}generated_files/png/{self.len_dir}/{self.file}_ss.png"

    @property
    def cm_status_json(self) -> str:
        return f"{self.arm_prefix}generated_files/cm/{self.len_dir}/{self.file}.cm.status.json"


def read_manifest(path: str | Path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def write_manifest(path: str | Path, entries: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.writelines(f"{entry}\n" for entry in entries)


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def run_checked(
    cmd: list[str],
    stdin_path: str | Path | None = None,
    stdout_path: str | Path | None = None,
    stderr_path: str | Path | None = None,
    cwd: str | Path | None = None,
) -> None:
    with contextlib.ExitStack() as stack:
        # Never inherit the caller's stdin. Some external tools (notably
        # refold.pl, whose <> falls back to STDIN once @ARGV is exhausted) will
        # otherwise block forever reading the user's terminal, stalling every
        # scheduler slot instead of failing.
        stdin_handle = (
            stack.enter_context(open(stdin_path, encoding="utf-8")) if stdin_path else subprocess.DEVNULL
        )
        stdout_handle = stack.enter_context(open(stdout_path, "w", encoding="utf-8")) if stdout_path else None
        stderr_handle = stack.enter_context(open(stderr_path, "w", encoding="utf-8")) if stderr_path else None
        subprocess.run(
            cmd,
            stdin=stdin_handle,
            stdout=stdout_handle,
            stderr=stderr_handle,
            cwd=cwd,
            text=True,
            check=True,
        )


def normalize_rnaalifold_side_output(outdir: Path, canonical_path: str | Path, suffix: str) -> None:
    canonical = Path(canonical_path)
    if canonical.exists():
        return

    prefix = canonical.name[: -len(suffix)]
    candidates = sorted(outdir.glob(f"{prefix}_*{suffix}"))
    if len(candidates) == 1:
        candidates[0].rename(canonical)
        return
    if len(candidates) > 1:
        raise FileExistsError(
            f"Multiple RNAalifold outputs matched {canonical.name}: "
            + ", ".join(candidate.name for candidate in candidates)
        )
    raise FileNotFoundError(
        f"Expected RNAalifold output not found and no renamed variant matched: {canonical.name}"
    )


def split_file_basenames_from_manifest(path: str | Path) -> list[str]:
    return [file[:-4] for file in read_manifest(path)]


def initial_alignment_input(input_alignment: str | None) -> list[str]:
    if not input_alignment:
        raise ValueError(
            "Missing required config value 'input_alignment'. "
            "Use RNAcs --input-alignment /path/to/input_alignment.{stk,aln}"
        )
    return [str(input_alignment)]


def initial_alignment_format_code(input_alignment: str | None) -> str:
    if not input_alignment:
        raise ValueError(
            "Missing required config value 'input_alignment'. "
            "Use RNAcs --input-alignment /path/to/input_alignment.{stk,aln}"
        )
    suffix = Path(input_alignment).suffix.lower()
    if suffix == ".aln":
        return "C"
    return "S"


def candidate_outputs_for_manifest(
    manifest_path: str | Path,
    wlen: str | int,
    path_getter: Callable[[CandidatePaths], str | list[str]],
    arm_prefix: str = "",
) -> list[str]:
    outputs: list[str] = []
    for file in split_file_basenames_from_manifest(manifest_path):
        value = path_getter(CandidatePaths(wlen=wlen, file=file, arm_prefix=arm_prefix))
        if isinstance(value, list):
            outputs.extend(value)
        else:
            outputs.append(value)
    return outputs
