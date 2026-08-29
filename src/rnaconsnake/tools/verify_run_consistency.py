from __future__ import annotations

import argparse
import hashlib
import json
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from rnaconsnake.workflow_helpers import analysis_root, read_manifest


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_window_lengths(run_dir: Path) -> list[int]:
    windows: set[int] = set()
    for base in (run_dir / "Lalifold", run_dir / "generated_files" / "stk"):
        if not base.is_dir():
            continue
        for child in base.iterdir():
            if child.is_dir() and child.name.startswith("len_"):
                try:
                    windows.add(int(child.name.split("_", 1)[1]))
                except ValueError:
                    continue
    return sorted(windows)


@dataclass(frozen=True)
class CandidateSnapshot:
    multistk_sha256: str
    split_manifest: list[str]
    split_hashes: dict[str, str]
    stk_manifest: list[str]
    stk_hashes: dict[str, str]


def build_snapshot(run_dir: Path, wlen: int) -> CandidateSnapshot:
    len_dir = f"len_{wlen}"
    multistk = run_dir / "Lalifold" / len_dir / f"RC_{wlen}_0001.stk"
    split_manifest_path = run_dir / "Lalifold" / len_dir / "split" / "manifest.txt"
    stk_manifest_path = run_dir / "generated_files" / "stk" / len_dir / "manifest.txt"

    if not multistk.is_file():
        raise FileNotFoundError(f"Missing RNALalifold multi-record output: {multistk}")
    if not split_manifest_path.is_file():
        raise FileNotFoundError(f"Missing split manifest: {split_manifest_path}")
    if not stk_manifest_path.is_file():
        raise FileNotFoundError(f"Missing cleaned alignment manifest: {stk_manifest_path}")

    split_manifest = read_manifest(split_manifest_path)
    stk_manifest = read_manifest(stk_manifest_path)

    split_hashes = {
        name: sha256_file(run_dir / "Lalifold" / len_dir / "split" / name) for name in split_manifest
    }
    stk_hashes = {
        name: sha256_file(run_dir / "generated_files" / "stk" / len_dir / name) for name in stk_manifest
    }

    return CandidateSnapshot(
        multistk_sha256=sha256_file(multistk),
        split_manifest=split_manifest,
        split_hashes=split_hashes,
        stk_manifest=stk_manifest,
        stk_hashes=stk_hashes,
    )


def diff_named_hashes(left: dict[str, str], right: dict[str, str], label: str) -> list[str]:
    messages: list[str] = []
    left_keys = set(left)
    right_keys = set(right)

    only_left = sorted(left_keys - right_keys)
    only_right = sorted(right_keys - left_keys)
    if only_left:
        messages.append(f"{label}: only in left: {', '.join(only_left[:10])}")
    if only_right:
        messages.append(f"{label}: only in right: {', '.join(only_right[:10])}")

    changed = sorted(name for name in left_keys & right_keys if left[name] != right[name])
    if changed:
        preview = ", ".join(changed[:10])
        suffix = " ..." if len(changed) > 10 else ""
        messages.append(f"{label}: content differs for {preview}{suffix}")
    return messages


def compare_snapshots(left: CandidateSnapshot, right: CandidateSnapshot, wlen: int) -> list[str]:
    messages: list[str] = []
    prefix = f"len_{wlen}"

    if left.multistk_sha256 != right.multistk_sha256:
        messages.append(f"{prefix}: RNALalifold multi-record output differs")
    if left.split_manifest != right.split_manifest:
        messages.append(f"{prefix}: split manifest differs")
    messages.extend(diff_named_hashes(left.split_hashes, right.split_hashes, f"{prefix}: split candidates"))
    if left.stk_manifest != right.stk_manifest:
        messages.append(f"{prefix}: cleaned alignment manifest differs")
    messages.extend(diff_named_hashes(left.stk_hashes, right.stk_hashes, f"{prefix}: cleaned alignments"))
    return messages


def comparison_payload(left_run: Path, right_run: Path, window_lengths: Iterable[int]) -> dict[str, object]:
    window_lengths = list(window_lengths)
    left_root = analysis_root(left_run)
    right_root = analysis_root(right_run)
    payload: dict[str, object] = {
        "left_run": str(left_run),
        "right_run": str(right_run),
        # Which layout was actually compared: the run directory, or the real
        # arm inside it. Recorded so a reader can tell the two apart.
        "analysis_root": {"left": str(left_root), "right": str(right_root)},
        "window_lengths": window_lengths,
        "identical": True,
        "differences": [],
    }
    differences: list[str] = []
    for wlen in window_lengths:
        differences.extend(
            compare_snapshots(build_snapshot(left_root, wlen), build_snapshot(right_root, wlen), wlen)
        )
    payload["differences"] = differences
    payload["identical"] = not differences
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare the deterministic candidate-generation outputs from two RNAConSnake runs. "
            "This intentionally ignores non-deterministic downstream metrics such as alifoldz z-scores."
        )
    )
    parser.add_argument("left_run", help="First run directory")
    parser.add_argument("right_run", help="Second run directory")
    parser.add_argument(
        "--maxbpspan",
        type=int,
        action="append",
        default=[],
        help="Restrict comparison to one or more window lengths. Repeat as needed.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a machine-readable JSON comparison result.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    left_run = Path(args.left_run).resolve()
    right_run = Path(args.right_run).resolve()

    # A calibrated run keeps every output under arms/real/. Discovering window
    # lengths in the run directory itself would find none there, compare
    # nothing, and report the two runs as identical.
    left_root = analysis_root(left_run)
    right_root = analysis_root(right_run)

    left_windows = discover_window_lengths(left_root)
    right_windows = discover_window_lengths(right_root)

    if args.maxbpspan:
        window_lengths = sorted(set(args.maxbpspan))
    else:
        if left_windows != right_windows:
            payload = {
                "left_run": str(left_run),
                "right_run": str(right_run),
                "identical": False,
                "differences": [
                    f"window lengths differ: left={left_windows}, right={right_windows}",
                ],
            }
            if args.json:
                print(json.dumps(payload, indent=2, sort_keys=True))
            else:
                for message in payload["differences"]:
                    print(message)
            return 1
        window_lengths = left_windows

    if not window_lengths:
        # Having compared nothing is not the same as having found no
        # difference, and must never be reported as agreement.
        payload = {
            "left_run": str(left_run),
            "right_run": str(right_run),
            "identical": False,
            "differences": [
                "no window lengths to compare: found no Lalifold/len_* or "
                f"generated_files/stk/len_* under {left_root} or {right_root}. "
                "Point this at completed run directories.",
            ],
        }
        if args.json:
            print(json.dumps(payload, indent=2, sort_keys=True))
        else:
            for message in payload["differences"]:
                print(message)
        return 1

    payload = comparison_payload(left_run, right_run, window_lengths)
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        if payload["identical"]:
            print("Deterministic candidate-generation outputs are identical.")
            print(f"Compared window lengths: {', '.join(str(wlen) for wlen in window_lengths)}")
            print(
                "Compared artifacts: RNALalifold multi-record outputs, split candidate "
                "manifests/files, cleaned alignment manifests/files."
            )
        else:
            for message in payload["differences"]:
                print(message)
    return 0 if payload["identical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
