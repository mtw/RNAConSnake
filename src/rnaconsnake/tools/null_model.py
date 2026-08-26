#!/usr/bin/env python3

"""Null-model alignment generation for the RNAConSnake calibration arm.

Two backends are wrapped behind one interface so that adding a third method
touches only this file:

``sissiz``
    SISSIz in alignment-simulation mode (``-s``).  Preserves the tree,
    dinucleotide content and the gap pattern.  This is the method the
    manuscript reports (Gesell & Washietl, BMC Bioinformatics 2008).
    SISSIz has no seed option and draws its own from the clock, so
    ``null.seed`` cannot reproduce a pool; pin the pool instead
    (:func:`adopt_pool`).
``rnazRandomizeAln``
    Column shuffling via the RNAz helper script.  Faster and weaker; kept as a
    smoke-test option and for CI.
``none``
    No null arm.  Never reaches this module: the workflow drops the whole
    ``arms/`` prefix instead.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from rnaconsnake.tools.alignment_io import (
    Alignment,
    alignment_diagnostics,
    gap_mask,
    mean_pairwise_identity,
    parse_clustal_blocks,
    read_stockholm_alignment,
    uses_rna_alphabet,
    write_clustal,
    write_stockholm_alignment,
)
from rnaconsnake.tools.stockholm_utils import parse_stockholm_records
from rnaconsnake.workflow_helpers import (
    REAL_ARM,
    arm_seed,
    perl_seed_env,
    perl_seeded_command,
    write_json,
)


MPI_DRIFT_WARNING = 0.05

# Seeds Perl's RNG before handing control to the RNAz helper script, so that a
# rerun with the same arm seed reproduces the same shuffled alignment.
PERL_SEED_BOOTSTRAP = "srand(shift); my $script = shift; do $script; die $@ if $@;"


class NullModelError(RuntimeError):
    pass


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run(cmd: list[str], stdout_path: str | Path | None = None, env: dict[str, str] | None = None) -> str:
    result = subprocess.run(
        cmd,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    if result.returncode != 0:
        raise NullModelError(
            f"null-model backend failed (exit {result.returncode}): {shlex.join(cmd)}\n"
            + (result.stderr or "").strip()
        )
    if stdout_path is not None:
        Path(stdout_path).write_text(result.stdout, encoding="utf-8")
    return result.stdout


def _perl_seeded_command(tokens: list[str], seed: int, args: list[str]) -> tuple[list[str], bool]:
    try:
        return perl_seeded_command(tokens, seed, args)
    except ValueError as error:
        raise NullModelError(
            "Empty command configured for the rnazRandomizeAln backend"
        ) from error


def _perl_env() -> dict[str, str]:
    return perl_seed_env()


def _align_names(source: Alignment, replicate: Alignment) -> tuple[Alignment, bool]:
    """Restore the source sequence names when the backend rewrote them."""
    if replicate.order == source.order:
        return replicate, False
    if len(replicate.order) != len(source.order):
        raise NullModelError(
            f"Simulated alignment has {len(replicate.order)} sequences, source has {len(source.order)}"
        )
    remapped = {
        new: replicate.seqs[old] for new, old in zip(source.order, replicate.order)
    }
    return Alignment(order=list(source.order), seqs=remapped), True


def validate_replicate(source: Alignment, replicate: Alignment, method: str) -> list[str]:
    warnings: list[str] = []
    if replicate.length != source.length:
        raise NullModelError(
            f"Simulated alignment length {replicate.length} != source length {source.length}"
        )

    source_masks = gap_mask(source)
    replicate_masks = gap_mask(replicate)
    if method == "sissiz":
        if source_masks != replicate_masks:
            raise NullModelError("SISSIz replicate did not preserve the source gap pattern")
    else:
        source_counts = [mask.count("-") for mask in source_masks]
        replicate_counts = [mask.count("-") for mask in replicate_masks]
        if source_counts != replicate_counts:
            raise NullModelError(
                "Shuffled replicate did not preserve per-sequence gap counts "
                f"({replicate_counts} != {source_counts})"
            )

    drift = abs(mean_pairwise_identity(replicate) - mean_pairwise_identity(source))
    if drift > MPI_DRIFT_WARNING:
        warnings.append(
            f"mean pairwise identity drifted by {drift:.3f} (> {MPI_DRIFT_WARNING})"
        )
    return warnings


def _simulate_sissiz(
    source: Alignment,
    replicates: int,
    command: list[str],
    workdir: Path,
) -> tuple[list[Alignment], dict[str, object]]:
    clustal_in = workdir / "source.aln"
    write_clustal(source, clustal_in)
    cmd = [*command, "-s", "-n", str(replicates), "--clustal"]
    if uses_rna_alphabet(source):
        cmd.append("--rna")
    cmd.append(str(clustal_in))
    stdout = _run(cmd)
    (workdir / "sissiz.out").write_text(stdout, encoding="utf-8")
    blocks = parse_clustal_blocks(stdout)
    if len(blocks) < replicates:
        raise NullModelError(
            f"SISSIz returned {len(blocks)} simulated alignments, expected {replicates}"
        )
    # SISSIz exposes no seed option and draws its own from the clock (two runs
    # inside the same second agree; a second apart they do not). The arm seed
    # therefore cannot reproduce a SISSIz pool -- reproducibility comes from
    # retaining or pinning the generated pool, whose SHA-256 is recorded.
    return blocks[:replicates], {
        "determinism": "unseeded",
        "seeded": False,
        "reproducibility": "retain_or_pin_pool",
    }


def _simulate_randomize(
    source: Alignment,
    replicates: int,
    command: list[str],
    workdir: Path,
    base_seed: int,
) -> tuple[list[Alignment], dict[str, object]]:
    clustal_in = workdir / "source.aln"
    write_clustal(source, clustal_in)
    blocks: list[Alignment] = []
    seeded_all = True
    seeds: list[int] = []
    for index in range(replicates):
        seed = arm_seed(f"null_{index:03d}", base_seed)
        assert seed is not None
        seeds.append(seed)
        cmd, seeded = _perl_seeded_command(command, seed, [str(clustal_in)])
        seeded_all = seeded_all and seeded
        stdout = _run(cmd, env=_perl_env())
        parsed = parse_clustal_blocks(stdout)
        if not parsed:
            raise NullModelError(
                f"rnazRandomizeAln produced no alignment for replicate {index}"
            )
        blocks.append(parsed[0])
    return blocks, {
        "determinism": "seeded" if seeded_all else "unseeded",
        "seeded": seeded_all,
        "replicate_seeds": seeds,
    }


def simulate_pool(
    source_path: str | Path,
    output_path: str | Path,
    metadata_path: str | Path,
    method: str,
    replicates: int,
    seed: int,
    sissiz_command: list[str],
    randomize_command: list[str],
    workdir: str | Path | None = None,
) -> dict[str, object]:
    """Generate ``replicates`` null alignments and write them as one Stockholm file."""
    if replicates < 1:
        raise NullModelError("simulate_pool requires at least one replicate")

    source = read_stockholm_alignment(source_path)
    if source.n_seq < 2:
        raise NullModelError(
            f"Null simulation needs at least 2 sequences, source has {source.n_seq}"
        )

    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(workdir) if workdir else Path(tmp)
        scratch.mkdir(parents=True, exist_ok=True)
        if method == "sissiz":
            blocks, backend_meta = _simulate_sissiz(source, replicates, sissiz_command, scratch)
        elif method == "rnazRandomizeAln":
            blocks, backend_meta = _simulate_randomize(
                source, replicates, randomize_command, scratch, seed
            )
        else:
            raise NullModelError(f"Unsupported null.method for simulation: {method!r}")

    warnings: list[str] = []
    if not backend_meta.get("seeded", False):
        warnings.append(
            f"the {method} backend is not seedable, so null.seed does not reproduce this "
            "pool; pin null_pool/pool.stk (see null.pool_file) to reproduce a calibration"
        )
    replicate_diagnostics: list[dict[str, object]] = []
    normalized: list[Alignment] = []
    for index, block in enumerate(blocks):
        aligned, remapped = _align_names(source, block)
        replicate_warnings = validate_replicate(source, aligned, method)
        arm = f"null_{index:03d}"
        for message in replicate_warnings:
            warnings.append(f"{arm}: {message}")
        diagnostics = alignment_diagnostics(aligned)
        diagnostics["arm"] = arm
        diagnostics["name_remapped"] = remapped
        replicate_diagnostics.append(diagnostics)
        normalized.append(aligned)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        for index, alignment in enumerate(normalized):
            handle.write("# STOCKHOLM 1.0\n")
            handle.write(f"#=GF ID null_{index:03d}\n")
            for name in alignment.order:
                handle.write(f"{name} {alignment.seqs[name]}\n")
            handle.write("//\n")

    metadata = {
        "method": method,
        "replicates": replicates,
        "seed": seed,
        "source": str(source_path),
        "source_diagnostics": alignment_diagnostics(source),
        "replicate_diagnostics": replicate_diagnostics,
        "pool_sha256": sha256_file(output_path),
        "warnings": warnings,
        **backend_meta,
    }
    write_json(metadata_path, metadata)
    for message in warnings:
        print(f"WARNING: null pool: {message}", file=sys.stderr)
    return metadata


def adopt_pool(
    pool_path: str | Path,
    source_path: str | Path,
    output_path: str | Path,
    metadata_path: str | Path,
    replicates: int,
    method: str,
) -> dict[str, object]:
    """Reuse a previously generated pool instead of simulating a new one.

    SISSIz cannot be seeded, so this is how a calibration is made reproducible:
    keep the pool that produced a result and pin it on the rerun. The pinned
    pool is validated against the current source alignment exactly as a freshly
    simulated one would be, so a pool from a different alignment is rejected
    rather than silently miscalibrating.
    """
    source = read_stockholm_alignment(source_path)
    records = parse_stockholm_records(pool_path)
    if len(records) < replicates:
        raise NullModelError(
            f"Pinned pool {pool_path} holds {len(records)} replicates, need {replicates}"
        )

    warnings: list[str] = []
    diagnostics: list[dict[str, object]] = []
    for index, record in enumerate(records[:replicates]):
        replicate = Alignment(order=list(record.seq_order), seqs=dict(record.seqs))
        aligned, remapped = _align_names(source, replicate)
        for message in validate_replicate(source, aligned, method):
            warnings.append(f"null_{index:03d}: {message}")
        entry = alignment_diagnostics(aligned)
        entry["arm"] = f"null_{index:03d}"
        entry["name_remapped"] = remapped
        diagnostics.append(entry)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(pool_path, output_path)
    metadata = {
        "method": method,
        "replicates": replicates,
        "source": str(source_path),
        "source_diagnostics": alignment_diagnostics(source),
        "replicate_diagnostics": diagnostics,
        "pool_sha256": sha256_file(output_path),
        "warnings": warnings,
        "determinism": "pinned",
        "seeded": False,
        "reproducibility": "pinned_pool",
        "pinned_from": str(pool_path),
    }
    write_json(metadata_path, metadata)
    return metadata


def make_arm_alignment(
    arm: str,
    source_path: str | Path,
    output_path: str | Path,
    pool_path: str | Path | None = None,
) -> None:
    """Populate ``arms/<arm>/input_alignment.stk``.

    The real arm gets a *copy* rather than a symlink: downstream tools
    occasionally rewrite alignments in place and a symlink would corrupt the
    source alignment.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    if arm == REAL_ARM:
        source = Path(source_path)
        if source.suffix.lower() == ".aln":
            blocks = parse_clustal_blocks(source.read_text(encoding="utf-8"))
            if len(blocks) != 1:
                raise NullModelError(
                    f"Expected exactly one Clustal alignment in {source}, found {len(blocks)}"
                )
            write_stockholm_alignment(blocks[0], output, identifier=REAL_ARM)
        else:
            shutil.copyfile(source, output)
        return

    if pool_path is None:
        raise NullModelError(f"Arm {arm!r} requires a simulated null pool")

    records = parse_stockholm_records(pool_path)
    wanted = [
        record
        for record in records
        if any(line.strip() == f"#=GF ID {arm}" for line in record.gf_lines)
    ]
    if len(wanted) != 1:
        raise NullModelError(
            f"Expected exactly one record with ID {arm} in {pool_path}, found {len(wanted)}"
        )
    record = wanted[0]
    with open(output, "w", encoding="utf-8") as handle:
        handle.write("# STOCKHOLM 1.0\n")
        handle.write(f"#=GF ID {arm}\n")
        for name in record.seq_order:
            handle.write(f"{name} {record.seqs[name]}\n")
        handle.write("//\n")


def cmd_simulate_pool(args: argparse.Namespace) -> int:
    if args.pool_file:
        adopt_pool(
            pool_path=args.pool_file,
            source_path=args.source,
            output_path=args.output,
            metadata_path=args.metadata,
            replicates=args.replicates,
            method=args.method,
        )
        return 0
    simulate_pool(
        source_path=args.source,
        output_path=args.output,
        metadata_path=args.metadata,
        method=args.method,
        replicates=args.replicates,
        seed=args.seed,
        sissiz_command=shlex.split(args.sissiz_command),
        randomize_command=shlex.split(args.randomize_command),
        workdir=args.workdir,
    )
    return 0


def cmd_make_arm(args: argparse.Namespace) -> int:
    make_arm_alignment(
        arm=args.arm,
        source_path=args.source,
        output_path=args.output,
        pool_path=args.pool,
    )
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="RNAConSnake null-model alignment generation.")
    sub = parser.add_subparsers(dest="command", required=True)

    pool = sub.add_parser("simulate-pool", help="Simulate all null replicates in one go.")
    pool.add_argument("--source", required=True)
    pool.add_argument("--output", required=True)
    pool.add_argument("--metadata", required=True)
    pool.add_argument("--method", required=True, choices=["sissiz", "rnazRandomizeAln"])
    pool.add_argument("--replicates", type=int, required=True)
    pool.add_argument("--seed", type=int, required=True)
    pool.add_argument("--sissiz-command", default="SISSIz")
    pool.add_argument("--randomize-command", default="rnazRandomizeAln.pl")
    pool.add_argument("--workdir")
    pool.add_argument(
        "--pool-file",
        help="Reuse this previously generated pool instead of simulating a new one.",
    )
    pool.set_defaults(func=cmd_simulate_pool)

    arm = sub.add_parser("make-arm", help="Materialise one arm's input alignment.")
    arm.add_argument("--arm", required=True)
    arm.add_argument("--source", required=True)
    arm.add_argument("--output", required=True)
    arm.add_argument("--pool")
    arm.set_defaults(func=cmd_make_arm)

    return parser


def main() -> int:
    args = build_parser().parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
