from __future__ import annotations

import contextlib
import json
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable


SUMMARY_FIELDS = [
    "wbn",
    "nrseq",
    "alilen",
    "maxcovarval",
    "maxcovarcount",
    "rscape_covary_count",
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
    def from_config(cls, config: dict[str, Any]) -> "WorkflowSettings":
        return cls(
            input_alignment=config.get("input_alignment"),
            maxbpspan=[int(value) for value in config.get("maxbpspan", [100, 200])],
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


@dataclass(frozen=True)
class CandidatePaths:
    wlen: str | int
    file: str

    @property
    def len_dir(self) -> str:
        return f"len_{self.wlen}"

    @property
    def split(self) -> str:
        return f"Lalifold/{self.len_dir}/split/{self.file}.stk"

    @property
    def orig(self) -> str:
        return f"generated_files/orig/{self.len_dir}/{self.file}.orig.stk"

    @property
    def remgap(self) -> str:
        return f"generated_files/remgap/{self.len_dir}/{self.file}_remgap.stk"

    @property
    def strip(self) -> str:
        return f"generated_files/strip/{self.len_dir}/{self.file}_stripped.stk"

    @property
    def stk(self) -> str:
        return f"generated_files/stk/{self.len_dir}/{self.file}.stk"

    @property
    def aln(self) -> str:
        return f"generated_files/aln/{self.len_dir}/{self.file}.aln"

    @property
    def rnaz_txt(self) -> str:
        return f"generated_files/rnaz/{self.len_dir}/{self.file}.rnaz.txt"

    @property
    def rnaz_json(self) -> str:
        return f"generated_files/rnaz/{self.len_dir}/{self.file}.rnaz.json"

    @property
    def alifoldz_txt(self) -> str:
        return f"generated_files/alifoldz/{self.len_dir}/{self.file}.alifoldz.txt"

    @property
    def alifoldz_json(self) -> str:
        return f"generated_files/alifoldz/{self.len_dir}/{self.file}.alifoldz.json"

    @property
    def rnaalifold_dir(self) -> str:
        return f"generated_files/rnaalifold/{self.len_dir}/{self.file}"

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
        return f"generated_files/refold/{self.len_dir}/{self.file}_refold.out"

    @property
    def refold_json(self) -> str:
        return f"generated_files/refold/{self.len_dir}/{self.file}.refold.json"

    @property
    def maxcovar_log(self) -> str:
        return f"generated_files/maxcovar/{self.len_dir}/{self.file}_alifoldmaxcovar.log"

    @property
    def maxcovar_json(self) -> str:
        return f"generated_files/maxcovar/{self.len_dir}/{self.file}.maxcovar.json"

    @property
    def rscape_power(self) -> str:
        return f"generated_files/rscape/{self.len_dir}/{self.file}.power"

    @property
    def rscape_json(self) -> str:
        return f"generated_files/rscape/{self.len_dir}/{self.file}.rscape.json"

    @property
    def rscape_sto_pdf(self) -> str:
        return f"generated_files/rscape/{self.len_dir}/{self.file}.sto.pdf"

    @property
    def summary_json(self) -> str:
        return f"generated_files/summary/{self.len_dir}/{self.file}.summary.json"

    @property
    def png_aln(self) -> str:
        return f"generated_files/png/{self.len_dir}/{self.file}_aln.png"

    @property
    def png_ss(self) -> str:
        return f"generated_files/png/{self.len_dir}/{self.file}_ss.png"

    @property
    def cm_status_json(self) -> str:
        return f"generated_files/cm/{self.len_dir}/{self.file}.cm.status.json"


def candidate_paths(wlen: str | int, file: str) -> CandidatePaths:
    return CandidatePaths(wlen=wlen, file=file)


def read_manifest(path: str | Path) -> list[str]:
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def write_manifest(path: str | Path, entries: list[str]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(f"{entry}\n")


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
        stdin_handle = stack.enter_context(open(stdin_path, encoding="utf-8")) if stdin_path else None
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


def split_file_basenames_from_manifest(path: str | Path) -> list[str]:
    return [file[:-4] for file in read_manifest(path)]


def initial_alignment_input(input_alignment: str | None) -> list[str]:
    if not input_alignment:
        raise ValueError(
            "Missing required config value 'input_alignment'. "
            "Use rnaconsnake-run --input-alignment /path/to/input_alignment.{stk,aln}"
        )
    return [str(input_alignment)]


def initial_alignment_format_code(input_alignment: str | None) -> str:
    if not input_alignment:
        raise ValueError(
            "Missing required config value 'input_alignment'. "
            "Use rnaconsnake-run --input-alignment /path/to/input_alignment.{stk,aln}"
        )
    suffix = Path(input_alignment).suffix.lower()
    if suffix == ".aln":
        return "C"
    return "S"


def candidate_outputs_for_manifest(
    manifest_path: str | Path,
    wlen: str | int,
    path_getter: Callable[[CandidatePaths], str | list[str]],
) -> list[str]:
    outputs: list[str] = []
    for file in split_file_basenames_from_manifest(manifest_path):
        value = path_getter(candidate_paths(wlen, file))
        if isinstance(value, list):
            outputs.extend(value)
        else:
            outputs.append(value)
    return outputs
