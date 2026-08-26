import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from rnaconsnake.workflow_helpers import (
    ARM_WILDCARD_PATTERN,
    REAL_ARM,
    CalibrationSettings,
    NullSettings,
    WorkflowSettings,
    CandidatePaths,
    arm_prefix_for,
    candidate_outputs_for_manifest,
    initial_alignment_format_code,
    initial_alignment_input as required_initial_alignment_input,
    derived_seed,
    normalize_rnaalifold_side_output,
    perl_seed_env,
    perl_seeded_command,
    read_json,
    read_manifest,
    run_checked,
    split_file_basenames_from_manifest,
    write_json,
    write_manifest,
)
from snakemake.io import expand


configfile: "config.yaml"

SETTINGS = WorkflowSettings.from_config(config)
NULL = NullSettings.from_config(config)
CALIBRATION = CalibrationSettings.from_config(config)
INPUT_ALIGNMENT = SETTINGS.input_alignment
INPUT_ALIGNMENT_FORMAT = initial_alignment_format_code(INPUT_ALIGNMENT)
MAXBPSPAN = SETTINGS.maxbpspan
LALIFOLD_THREADS = SETTINGS.lalifold_threads
REMOVE_GAPONLY_GAPRATIO = SETTINGS.remove_gaponly_gapratio
REMOVE_GAPONLY_MAX_N = SETTINGS.remove_gaponly_max_n
DO_CM = SETTINGS.do_cm
DO_LOCARNATE = SETTINGS.do_locarnate
DO_PNG = SETTINGS.do_png
DO_RSCAPE = SETTINGS.do_rscape
RNAZ_NO_SHUFFLE = SETTINGS.rnaz_no_shuffle
CM_RNAZ_THRESHOLD = SETTINGS.cm_rnaz_threshold
CM_ALIFOLDZ_THRESHOLD = SETTINGS.cm_alifoldz_threshold

# --- null-model calibration arm ------------------------------------------
#
# Every rule downstream of alignment generation is traversed by *both* arms.
# There are no parallel rules: the only difference between the real arm and a
# null replicate is the ``arms/{arm}/`` path prefix, which collapses to the
# empty string when the null arm is disabled.  That is what makes
# ``null.method: none`` reproduce the pre-calibration outputs byte-for-byte,
# and what guarantees the two arms cannot silently diverge.
NULL_ARM_ENABLED = NULL.enabled
ARMS = NULL.arms()
ARM_PREFIX = "arms/{arm}/" if NULL_ARM_ENABLED else ""
NULL_POOL = "null_pool/pool.stk"
NULL_POOL_METADATA = "null_pool/pool.json"
CALIBRATION_DIR = "results/calibration"
BENCHMARK_DIR = "results/benchmark"
VERSIONS_FILE = "results/versions.yaml"
# AlifoldZ shuffles internally, so it is the expensive branch.  When two-stage
# mode is on it runs only on stage-one survivors -- in *both* arms, otherwise
# the calibration would be invalid -- and the resulting FDR is conditional on
# passing stage one.  That conditionality is recorded in summary.json and in
# the funnel header.
TWO_STAGE = NULL_ARM_ENABLED and NULL.two_stage
STAGE1_RNAZ_PROB = CALIBRATION.stage1_rnaz_prob
# alifoldz.pl shuffles internally and exposes no seed, so its z-scores vary
# between runs. Seeding Perl's RNG per candidate makes a run reproducible,
# which the calibration needs: q-values derive from these scores.
ALIFOLDZ_SEED = config.get("alifoldz_seed")
DEREPLICATE_METHOD = CALIBRATION.dereplicate_method
PAIR_CONTAINMENT = CALIBRATION.pair_containment
EMIT_VERSIONS = bool(config.get("emit_versions", False)) or NULL_ARM_ENABLED
BENCHMARK_TRUTH_SETTING = config.get(
    "benchmark_truth", "resources/benchmark/flavivirus_elements.tsv"
)
BENCHMARK_ALIGNMENT = config.get("benchmark_alignment")
BENCHMARK_MIN_OVERLAP = float(config.get("benchmark_min_overlap_fraction", 0.5))
BENCHMARK_ALLOW_UNCURATED = bool(config.get("benchmark_allow_uncurated", False))

wildcard_constraints:
    wlen=r"\d+",
    arm=ARM_WILDCARD_PATTERN


def resolve_benchmark_truth(setting):
    """Locate the curated truth file.

    The workflow usually runs with ``--directory <run_dir>``, so a
    repository-relative default would not resolve there. Look in the run
    directory first (a run-local override wins), then next to the snakefile,
    then in the installed package.
    """
    candidate = Path(setting)
    if candidate.is_absolute():
        return str(candidate)
    if candidate.exists():
        return str(candidate)

    basedir = Path(workflow.basedir)
    for root in [basedir, basedir.parent, basedir.parents[1] if len(basedir.parents) > 1 else basedir]:
        located = root / candidate
        if located.is_file():
            return str(located)

    packaged = basedir.parent / "resources" / "benchmark" / candidate.name
    if packaged.is_file():
        return str(packaged)
    # Fall through unchanged: the rule will fail with a missing-input error
    # naming the path the user configured.
    return str(candidate)


BENCHMARK_TRUTH = resolve_benchmark_truth(BENCHMARK_TRUTH_SETTING)


def A(path):
    """Prefix a workflow path with the current arm directory.

    Returns ``path`` unchanged when the null arm is disabled, so the legacy
    output layout is preserved exactly.
    """
    return ARM_PREFIX + path


def arm_expand(pattern, **kwargs):
    if NULL_ARM_ENABLED:
        return expand(pattern, arm=ARMS, **kwargs)
    return expand(pattern, **kwargs)


def arm_prefix_of(wildcards):
    return arm_prefix_for(getattr(wildcards, "arm", None) if NULL_ARM_ENABLED else None)


def checkpoint_wildcards(wildcards):
    keys = {"wlen": wildcards.wlen}
    if NULL_ARM_ENABLED:
        keys["arm"] = wildcards.arm
    return keys


def stage_one_pass(rnaz_json_path, threshold):
    """Cheap stage-one gate: did this candidate clear the RNAz screen?"""
    try:
        return float(read_json(rnaz_json_path).get("rnazprob", "-1")) >= threshold
    except (OSError, ValueError):
        return False


def paths_for(wildcards):
    return CandidatePaths(
        wlen=wildcards.wlen,
        file=wildcards.file,
        arm_prefix=arm_prefix_of(wildcards),
    )


def command_tokens(name, default):
    return SETTINGS.command_tokens(name, default)


def write_output_manifest(output_path, input_paths):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    write_manifest(output_path, [os.path.basename(p) for p in input_paths])


def split_file_basenames(wildcards):
    manifest = checkpoints.split_stockholm.get(**checkpoint_wildcards(wildcards)).output[1]
    return split_file_basenames_from_manifest(manifest)


def initial_alignment_input(wildcards):
    if NULL_ARM_ENABLED:
        return [f"arms/{wildcards.arm}/input_alignment.stk"]
    return required_initial_alignment_input(INPUT_ALIGNMENT)


def _candidate_outputs(wildcards, path_getter):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(**checkpoint_wildcards(wildcards)).output[1],
        wildcards.wlen,
        path_getter,
        arm_prefix=arm_prefix_of(wildcards),
    )


def orig_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.orig)


def remgap_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.remgap)


def strip_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.strip)


def stk_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.stk)


def cm_status_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.cm_status_json)


def summary_json_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.summary_json)


def rscape_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: paths.rscape_json)


def png_outputs(wildcards):
    return _candidate_outputs(wildcards, lambda paths: [paths.png_aln, paths.png_ss])


def arm_alignment_input(wildcards):
    if wildcards.arm == REAL_ARM:
        return required_initial_alignment_input(INPUT_ALIGNMENT)
    return [NULL_POOL]


def calibration_arm_inputs():
    """``(arm, wlen) -> per-arm summary CSV`` consumed by the calibration step."""
    return {
        (arm, wlen): f"arms/{arm}/generated_files/summary/len_{wlen}/RNAConSnake.log.csv"
        for arm in ARMS
        for wlen in MAXBPSPAN
    }


def pipeline_targets():
    targets = [
        arm_expand(A("Lalifold/len_{wlen}/RC_{wlen}_0001.stk"), wlen=MAXBPSPAN),
        arm_expand(A("Lalifold/len_{wlen}/split/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/orig/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/remgap/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/strip/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/stk/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/rscape/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/cm/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/summary/len_{wlen}/RNAConSnake.log"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/summary/len_{wlen}/RNAConSnake.log.csv"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/summary/len_{wlen}/RNAConSnake.md"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/summary/len_{wlen}/RNAConSnake.nr.csv"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/summary/len_{wlen}/RNAConSnake.nr.json"), wlen=MAXBPSPAN),
        arm_expand(A("generated_files/alignment_screenability.tsv")),
    ]
    if DO_PNG:
        targets.append(arm_expand(A("generated_files/png/len_{wlen}/manifest.txt"), wlen=MAXBPSPAN))
    if NULL_ARM_ENABLED:
        targets.append(
            [
                f"{CALIBRATION_DIR}/funnel.tsv",
                f"{CALIBRATION_DIR}/qvalues.tsv",
                f"{CALIBRATION_DIR}/score_dists.tsv",
                f"{CALIBRATION_DIR}/summary.json",
            ]
        )
    if EMIT_VERSIONS:
        targets.append([VERSIONS_FILE])
    return targets


rule all:
    input:
        pipeline_targets()


rule RNALalifold:
    input:
        initial_alignment_input
    output:
        stdout=A("Lalifold/len_{wlen}/RNALalifold.out"),
        stderr=A("Lalifold/len_{wlen}/RNALalifold.err"),
        multistk=A("Lalifold/len_{wlen}/RC_{wlen}_0001.stk")
    params:
        cmd=SETTINGS.tools.get("rnalalifold", "RNALalifold"),
        input_abs=lambda wildcards, input: os.path.abspath(input[0]),
        input_format=lambda wildcards: "S" if NULL_ARM_ENABLED else INPUT_ALIGNMENT_FORMAT,
        outdir=lambda wildcards, output: os.path.dirname(output.multistk)
    threads:
        LALIFOLD_THREADS
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        {params.cmd} \
            -L {wildcards.wlen} \
            --aln-stk \
            --id-prefix RC_{wildcards.wlen} \
            --cfactor 0.6 --nfactor 0.5 \
            -r \
            --csv \
            -f {params.input_format} \
            < {params.input_abs} > RNALalifold.out 2> RNALalifold.err
        """


checkpoint split_stockholm:
    input:
        A("Lalifold/len_{wlen}/RC_{wlen}_0001.stk")
    output:
        directory(A("Lalifold/len_{wlen}/split")),
        A("Lalifold/len_{wlen}/split/manifest.txt")
    log:
        out=A("Lalifold/len_{wlen}/split/split.out"),
        err=A("Lalifold/len_{wlen}/split/split.err")
    threads:
        1
    shell:
        """
        mkdir -p {output[0]}
        cd {output[0]}
        python3 -m rnaconsnake.tools.split_stockholm -a ../RC_{wildcards.wlen}_0001.stk > split.out 2> split.err
        find . -maxdepth 1 -type f -name 'RC_{wildcards.wlen}_*.stk' -print | sed 's#^\./##' | sort > manifest.txt
        """


rule preprocess_alignment_file:
    input:
        A("Lalifold/len_{wlen}/split/{file}.stk")
    output:
        orig=A("generated_files/orig/len_{wlen}/{file}.orig.stk"),
        remgap=A("generated_files/remgap/len_{wlen}/{file}_remgap.stk"),
        strip=A("generated_files/strip/len_{wlen}/{file}_stripped.stk"),
        stk=A("generated_files/stk/len_{wlen}/{file}.stk")
    threads:
        1
    run:
        paths = paths_for(wildcards)
        Path(output.orig).parent.mkdir(parents=True, exist_ok=True)
        Path(output.remgap).parent.mkdir(parents=True, exist_ok=True)
        Path(output.strip).parent.mkdir(parents=True, exist_ok=True)
        Path(output.stk).parent.mkdir(parents=True, exist_ok=True)

        run_checked(["cp", input[0], paths.orig])
        run_checked(
            command_tokens("remove_gaponly", "python3 -m rnaconsnake.tools.remove_gaponly")
            + ["-a", input[0], "-i", "stockholm", "-r", str(REMOVE_GAPONLY_GAPRATIO), "-n", str(REMOVE_GAPONLY_MAX_N)],
            stdout_path=paths.remgap,
            stderr_path=os.devnull,
        )
        run_checked(
            command_tokens("strip_aln", "python3 -m rnaconsnake.tools.strip_aln")
            + ["-a", paths.remgap, "-f", "S", "--nosingle"],
            stdout_path=paths.strip,
        )
        run_checked(["cp", paths.strip, paths.stk])


rule orig_manifest:
    input: orig_outputs
    output: A("generated_files/orig/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule remgap_manifest:
    input: remgap_outputs
    output: A("generated_files/remgap/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule strip_manifest:
    input: strip_outputs
    output: A("generated_files/strip/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule stk_manifest:
    input: stk_outputs
    output: A("generated_files/stk/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule analyze_alignment_file:
    input:
        A("generated_files/stk/len_{wlen}/{file}.stk")
    output:
        aln=A("generated_files/aln/len_{wlen}/{file}.aln"),
        rnaz_txt=A("generated_files/rnaz/len_{wlen}/{file}.rnaz.txt"),
        rnaz_metrics=A("generated_files/rnaz/len_{wlen}/{file}.rnaz.json"),
        alifoldz_txt=A("generated_files/alifoldz/len_{wlen}/{file}.alifoldz.txt"),
        alifoldz_metrics=A("generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json")
    run:
        paths = paths_for(wildcards)
        Path(paths.aln).parent.mkdir(parents=True, exist_ok=True)
        Path(paths.rnaz_txt).parent.mkdir(parents=True, exist_ok=True)
        Path(paths.alifoldz_txt).parent.mkdir(parents=True, exist_ok=True)

        run_checked(
            command_tokens("eslreformat", "esl-reformat") + ["clustal", input[0]],
            stdout_path=paths.aln,
        )

        cmd = command_tokens("rnaz", "RNAz") + ["-d"]
        if RNAZ_NO_SHUFFLE:
            cmd.append("-n")
        if DO_LOCARNATE:
            cmd.append("-l")
        cmd.append(paths.aln)
        run_checked(cmd, stdout_path=paths.rnaz_txt)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + ["extract-rnaz", "--input", paths.rnaz_txt, "--output", paths.rnaz_json]
        )

        # Two-stage design: AlifoldZ shuffles internally and is the expensive
        # branch, so it runs only on stage-one survivors.  The gate is applied
        # identically in the real and null arms; any divergence here would
        # invalidate the calibration.  Skipped candidates get a non-numeric
        # AlifoldZ score, which the calibration step counts as "did not pass".
        if TWO_STAGE and not stage_one_pass(paths.rnaz_json, STAGE1_RNAZ_PROB):
            Path(paths.alifoldz_txt).write_text(
                "# alifoldz skipped: candidate did not pass stage one "
                f"(RNAz class probability < {STAGE1_RNAZ_PROB})\n",
                encoding="utf-8",
            )
            write_json(paths.alifoldz_json, {"alifoldzscore": "NA"})
            return

        # alifoldz.pl shells out to `RNAalifold ... <tmpfile>` without --noPS
        # (alifoldz.pl:239), so RNAalifold drops an undeclared "alirna.ps" into
        # the *current* directory -- once per shuffle, from every candidate in
        # parallel, all racing on the same path in the run root. Give each
        # invocation a private scratch cwd so the stray plot is confined and
        # disposed of with it.
        alifoldz_cmd = command_tokens("alifoldz", "alifoldz.pl")
        alifoldz_env = None
        if ALIFOLDZ_SEED is not None:
            # A per-candidate seed: parallel candidates stay independent, and
            # each reproduces on a rerun.
            alifoldz_cmd, seeded = perl_seeded_command(
                alifoldz_cmd,
                derived_seed(int(ALIFOLDZ_SEED), wildcards.file),
                ["-f", "-t", "0.0"],
            )
            if seeded:
                alifoldz_env = perl_seed_env()
        else:
            alifoldz_cmd = alifoldz_cmd + ["-f", "-t", "0.0"]

        with tempfile.TemporaryDirectory(prefix="alifoldz-") as scratch:
            with open(paths.aln, encoding="utf-8") as stdin_handle:
                result = subprocess.run(
                    alifoldz_cmd,
                    stdin=stdin_handle,
                    cwd=scratch,
                    env=alifoldz_env,
                    capture_output=True,
                    text=True,
                    check=False,
                )

        with open(paths.alifoldz_txt, "w", encoding="utf-8") as handle:
            if result.stdout:
                handle.write(result.stdout)
            if result.returncode != 0:
                if result.stdout and not result.stdout.endswith("\n"):
                    handle.write("\n")
                stderr = (result.stderr or "").strip()
                handle.write(f"# alifoldz failed with exit code {result.returncode}\n")
                if stderr:
                    handle.write(f"# stderr: {stderr}\n")

        if result.returncode == 0:
            run_checked(
                command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
                + ["extract-alifoldz", "--input", paths.alifoldz_txt, "--output", paths.alifoldz_json]
            )
        else:
            # Not 0.0: a numeric-looking fallback is indistinguishable from a
            # real, unremarkable z-score and would silently enter the summary
            # tables and the FDR calibration. "NA" fails the AlifoldZ filter.
            write_json(paths.alifoldz_json, {"alifoldzscore": "NA"})


rule run_post_rnaalifold_file:
    input:
        A("generated_files/stk/len_{wlen}/{file}.stk")
    output:
        stdout=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.alifold.out"),
        stk=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"),
        ali_out=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ali.out"),
        dp_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_dp.ps"),
        aln_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.ps"),
        aln_eps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.eps"),
        aln_pdf=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.pdf"),
        ss_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.ps"),
        ss_eps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.eps"),
        ss_pdf=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.pdf")
    run:
        paths = paths_for(wildcards)
        outdir = Path(paths.rnaalifold_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        run_checked(
            command_tokens("rnaalifold", "RNAalifold")
            + [
                "-t4",
                "--aln",
                "--color",
                "-r",
                "--cfactor",
                "0.6",
                "--nfactor",
                "0.5",
                "-p",
                "--aln-EPS-cols=200",
                "--aln-stk",
                "-f",
                "S",
                "--id-prefix",
                wildcards.file,
            ],
            stdin_path=input[0],
            stdout_path=paths.rnaalifold_stdout,
            cwd=str(outdir),
        )
        # ViennaRNA can emit an undeclared default structure plot named
        # "alirna.ps". It is not part of RNAConSnake's output contract.
        # Only this job's own output directory is cleaned: reaching into the
        # run root from a per-candidate rule would race with the other
        # candidates still running.
        stray = outdir / "alirna.ps"
        if stray.exists():
            stray.unlink()
        default_stk = outdir / "RNAalifold_results.stk"
        if default_stk.exists() and not Path(paths.rnaalifold_stk).exists():
            default_stk.rename(paths.rnaalifold_stk)
        normalize_rnaalifold_side_output(outdir, paths.ali_out, "_ali.out")
        normalize_rnaalifold_side_output(outdir, paths.dp_ps, "_dp.ps")
        normalize_rnaalifold_side_output(outdir, paths.aln_ps, "_aln.ps")
        normalize_rnaalifold_side_output(outdir, paths.ss_ps, "_ss.ps")
        run_checked(command_tokens("ps2eps", "ps2eps") + [Path(paths.aln_ps).name], cwd=str(outdir))
        run_checked(command_tokens("epstopdf", "epstopdf") + [Path(paths.aln_eps).name], cwd=str(outdir))
        run_checked(command_tokens("ps2eps", "ps2eps") + [Path(paths.ss_ps).name], cwd=str(outdir))
        run_checked(command_tokens("epstopdf", "epstopdf") + [Path(paths.ss_eps).name], cwd=str(outdir))
        for required in [
            paths.rnaalifold_stk,
            paths.ali_out,
            paths.dp_ps,
            paths.aln_ps,
            paths.aln_eps,
            paths.aln_pdf,
            paths.ss_ps,
            paths.ss_eps,
            paths.ss_pdf,
        ]:
            if not Path(required).exists():
                raise FileNotFoundError(f"Expected RNAalifold output missing: {required}")


rule render_pngs_file:
    input:
        aln_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.ps"),
        ss_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.ps")
    output:
        aln_png=A("generated_files/png/len_{wlen}/{file}_aln.png"),
        ss_png=A("generated_files/png/len_{wlen}/{file}_ss.png")
    run:
        paths = paths_for(wildcards)
        os.makedirs(os.path.dirname(paths.png_aln), exist_ok=True)
        run_checked(command_tokens("magick", "magick") + [input.aln_ps, paths.png_aln])
        run_checked(command_tokens("magick", "magick") + [input.ss_ps, paths.png_ss])


rule png_manifest:
    input: png_outputs if DO_PNG else lambda wildcards: []
    output: A("generated_files/png/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule reformat_rnaalifold_results_file:
    input:
        A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk")
    output:
        A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln")
    run:
        run_checked(
            command_tokens("eslreformat", "esl-reformat") + ["clustal", input[0]],
            stdout_path=output[0],
        )


rule clean_rnaalifold_clustal_file:
    input:
        A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln")
    output:
        backup=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln~"),
        cleaned=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.cleaned.aln")
    run:
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "clean-clustal",
                "--input",
                input[0],
                "--backup",
                output.backup,
                "--output",
                output.cleaned,
            ]
        )


rule run_refold_file:
    """Refold each sequence under the consensus structure.

    One tool, not the ``refold.pl | RNAfold -C`` pipe it replaces: the
    constraints and the constrained fold both come from
    ``rnaconsnake.tools.refold``, which needs only the ViennaRNA Python
    bindings. It fails loudly on an alignment it cannot read, where the Perl
    script exited 0 having silently parsed nothing.
    """
    input:
        aln=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.cleaned.aln"),
        dp_ps=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_dp.ps")
    output:
        A("generated_files/refold/len_{wlen}/{file}_refold.out")
    threads:
        1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        run_checked(
            command_tokens("refold", "python3 -m rnaconsnake.tools.refold")
            + [
                "--alignment", input.aln,
                "--consensus", input.dp_ps,
                "--output", output[0],
            ]
        )


rule extract_refold_metrics_file:
    input:
        refold=A("generated_files/refold/len_{wlen}/{file}_refold.out"),
        stk=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk")
    output:
        A("generated_files/refold/len_{wlen}/{file}.refold.json")
    run:
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "extract-refold",
                "--rnaalifold-stk",
                input.stk,
                "--output",
                output[0],
            ]
        )


rule run_maxcovar_file:
    input:
        A("generated_files/rnaalifold/len_{wlen}/{file}/{file}_ali.out")
    output:
        log=A("generated_files/maxcovar/len_{wlen}/{file}_alifoldmaxcovar.log"),
        metrics=A("generated_files/maxcovar/len_{wlen}/{file}.maxcovar.json")
    run:
        paths = paths_for(wildcards)
        os.makedirs(os.path.dirname(paths.maxcovar_log), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "run-maxcovar",
                "--ali-out",
                input[0],
                "--log",
                paths.maxcovar_log,
                "--output",
                paths.maxcovar_json,
            ]
        )


rule run_rscape_file:
    input:
        A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk")
    output:
        power=A("generated_files/rscape/len_{wlen}/{file}.power"),
        metrics=A("generated_files/rscape/len_{wlen}/{file}.rscape.json")
    run:
        paths = paths_for(wildcards)
        os.makedirs(os.path.dirname(paths.rscape_power), exist_ok=True)
        if not DO_RSCAPE:
            Path(paths.rscape_power).write_text("# R-scape disabled\n", encoding="utf-8")
            write_json(paths.rscape_json, {"rscape_covary_count": "NA"})
        else:
            outdir = Path(paths.rscape_power).parent
            workdir = outdir / f".{wildcards.file}.rscape"
            workdir.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(
                command_tokens("rscape", "R-scape") + ["-s", os.path.abspath(input[0])],
                cwd=str(workdir),
                stdin=subprocess.DEVNULL,
                text=True,
                capture_output=True,
                check=False,
            )
            (workdir / "rscape.stdout").write_text(result.stdout or "", encoding="utf-8")
            (workdir / "rscape.stderr").write_text(result.stderr or "", encoding="utf-8")
            power_candidates = sorted(workdir.glob("*.power"))
            if len(power_candidates) == 1:
                power_candidates[0].replace(paths.rscape_power)
                run_checked(
                    command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
                    + [
                        "extract-rscape",
                        "--input",
                        paths.rscape_power,
                        "--output",
                        paths.rscape_json,
                    ]
                )
            elif len(power_candidates) == 0:
                stdout_text = (workdir / "rscape.stdout").read_text(encoding="utf-8") if (workdir / "rscape.stdout").exists() else ""
                stderr_text = (workdir / "rscape.stderr").read_text(encoding="utf-8") if (workdir / "rscape.stderr").exists() else ""
                Path(paths.rscape_power).write_text(
                    "# R-scape produced no .power output\n"
                    + f"# exit code: {result.returncode}\n"
                    + (stdout_text if stdout_text else "")
                    + ("\n# stderr\n" + stderr_text if stderr_text else ""),
                    encoding="utf-8",
                )
                write_json(paths.rscape_json, {"rscape_covary_count": "0" if "Number of covarying pairs = 0" in stdout_text else ""})
            else:
                candidates = ", ".join(path.name for path in power_candidates)
                raise FileNotFoundError(f"Could not uniquely identify R-scape .power output in {workdir}: {candidates}")
            sto_pdf = Path(paths.rscape_sto_pdf)
            if sto_pdf.exists():
                sto_pdf.unlink()
            metrics = read_json(paths.rscape_json)
            try:
                covary_count = int(metrics.get("rscape_covary_count", "") or 0)
            except ValueError:
                covary_count = 0
            if covary_count > 0:
                sto_candidates = sorted(workdir.glob("*.sto.pdf"))
                if len(sto_candidates) == 1:
                    shutil.copy2(sto_candidates[0], sto_pdf)
            else:
                shutil.rmtree(workdir, ignore_errors=True)


rule rscape_manifest:
    input: rscape_outputs
    output: A("generated_files/rscape/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule build_cm_file:
    input:
        rnaz=A("generated_files/rnaz/len_{wlen}/{file}.rnaz.json"),
        alifoldz=A("generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json"),
        stk=A("generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk")
    output:
        A("generated_files/cm/len_{wlen}/{file}.cm.status.json")
    run:
        paths = paths_for(wildcards)
        outdir = Path(paths.cm_status_json).parent
        outdir.mkdir(parents=True, exist_ok=True)

        rnaz = read_json(input.rnaz)
        alifoldz = read_json(input.alifoldz)

        try:
            rnazprob = float(rnaz.get("rnazprob", "-1"))
        except ValueError:
            rnazprob = -1.0
        try:
            alifoldzscore = float(alifoldz.get("alifoldzscore", "0"))
        except ValueError:
            alifoldzscore = 0.0

        status = {"built": False, "reason": "disabled"}
        build_reason = None
        if DO_CM:
            if rnazprob >= CM_RNAZ_THRESHOLD:
                build_reason = f"rnazprob={rnazprob}"
            elif alifoldzscore <= CM_ALIFOLDZ_THRESHOLD:
                build_reason = f"alifoldzscore={alifoldzscore}"

            if build_reason is not None:
                cm_base = str(outdir / paths.file)
                run_checked(
                    command_tokens("cmbuild", "cmbuild") + [f"{cm_base}.cm", input.stk],
                    stdout_path=f"{cm_base}.cmbuild.out",
                    stderr_path=f"{cm_base}.cmbuild.err",
                )
                run_checked(
                    command_tokens("cmcalibrate", "cmcalibrate") + [f"{cm_base}.cm"],
                    stdout_path=f"{cm_base}.cmcalibrate.out",
                    stderr_path=f"{cm_base}.cmcalibrate.err",
                )
                status = {"built": True, "reason": build_reason, "cm": f"{cm_base}.cm"}
            else:
                status = {"built": False, "reason": "threshold_not_met"}

        write_json(paths.cm_status_json, status)


rule cm_manifest:
    input: cm_status_outputs
    output: A("generated_files/cm/len_{wlen}/manifest.txt")
    run: write_output_manifest(output[0], input)


rule combine_summary_metrics_file:
    input:
        rnaz=A("generated_files/rnaz/len_{wlen}/{file}.rnaz.json"),
        alifoldz=A("generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json"),
        refold=A("generated_files/refold/len_{wlen}/{file}.refold.json"),
        maxcov=A("generated_files/maxcovar/len_{wlen}/{file}.maxcovar.json"),
        rscape=A("generated_files/rscape/len_{wlen}/{file}.rscape.json")
    output:
        A("generated_files/summary/len_{wlen}/{file}.summary.json")
    run:
        paths = paths_for(wildcards)
        os.makedirs(os.path.dirname(paths.summary_json), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "combine-summary",
                "--wbn",
                wildcards.file,
                "--output",
                paths.summary_json,
                input.rnaz,
                input.alifoldz,
                input.refold,
                input.maxcov,
                input.rscape,
            ]
        )


rule summary_logs:
    input:
        summary_json_outputs
    output:
        log=A("generated_files/summary/len_{wlen}/RNAConSnake.log"),
        csv=A("generated_files/summary/len_{wlen}/RNAConSnake.log.csv")
    run:
        os.makedirs(os.path.dirname(output.log), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "write-summary-outputs",
                "--label",
                f"len_{wildcards.wlen}",
                "--log",
                output.log,
                "--csv",
                output.csv,
                *sorted(input),
            ]
        )


# --- null-model calibration arm -------------------------------------------


if NULL_ARM_ENABLED:

    rule simulate_null_pool:
        """Simulate every null replicate in one backend invocation.

        SISSIz is deterministic for a given input alignment, so all replicates
        come from a single ``-n <replicates>`` run; asking for replicate *k*
        separately would cost O(k) simulations for the same result.
        """
        input:
            required_initial_alignment_input(INPUT_ALIGNMENT)
        output:
            pool=NULL_POOL,
            metadata=NULL_POOL_METADATA
        log:
            "logs/simulate_null_pool.log"
        params:
            method=NULL.method,
            replicates=NULL.replicates,
            seed=NULL.seed,
            pool_file=NULL.pool_file
        threads:
            1
        run:
            os.makedirs(os.path.dirname(output.pool), exist_ok=True)
            os.makedirs(os.path.dirname(log[0]), exist_ok=True)
            source = input[0]
            # The backends read Stockholm; convert a Clustal input first so
            # that both arms start from the same normalised alignment.
            if Path(source).suffix.lower() == ".aln":
                source = os.path.join(os.path.dirname(output.pool), "source.stk")
                run_checked(
                    command_tokens("null_model", "python3 -m rnaconsnake.tools.null_model")
                    + ["make-arm", "--arm", REAL_ARM, "--source", input[0], "--output", source]
                )
            run_checked(
                command_tokens("null_model", "python3 -m rnaconsnake.tools.null_model")
                + [
                    "simulate-pool",
                    "--source", source,
                    "--output", output.pool,
                    "--metadata", output.metadata,
                    "--method", params.method,
                    "--replicates", str(params.replicates),
                    "--seed", str(params.seed),
                    "--sissiz-command", SETTINGS.tools.get("sissiz", "SISSIz"),
                    "--randomize-command", SETTINGS.tools.get("rnaz_randomize_aln", "rnazRandomizeAln.pl"),
                    "--workdir", os.path.join(os.path.dirname(output.pool), "work"),
                ]
                + (["--pool-file", str(params.pool_file)] if params.pool_file else []),
                stderr_path=log[0],
            )

    rule make_arm_alignment:
        """Populate one arm's input alignment.

        This is the *only* point at which the arms differ.  The real arm gets a
        copy (never a symlink: downstream tools occasionally rewrite alignments
        in place, and a symlink would corrupt the source alignment); a null arm
        gets its replicate from the simulated pool.

        Per-arm seeds are derived by ``NullSettings.arm_seed`` and applied in
        ``simulate_null_pool``, where the randomness actually lives; this rule
        only selects an already-generated replicate.
        """
        input:
            arm_alignment_input
        output:
            "arms/{arm}/input_alignment.stk"
        log:
            "logs/{arm}/make_arm_alignment.log"
        threads:
            1
        run:
            os.makedirs(os.path.dirname(log[0]), exist_ok=True)
            cmd = command_tokens("null_model", "python3 -m rnaconsnake.tools.null_model") + [
                "make-arm",
                "--arm", wildcards.arm,
                "--source", str(INPUT_ALIGNMENT),
                "--output", output[0],
            ]
            if wildcards.arm != REAL_ARM:
                cmd += ["--pool", input[0]]
            run_checked(cmd, stderr_path=log[0])

    rule calibrate:
        """Aggregate all arms into empirical FDR / q-values and a filter funnel."""
        input:
            summaries=sorted(calibration_arm_inputs().values()),
            metadata=NULL_POOL_METADATA
        output:
            funnel=f"{CALIBRATION_DIR}/funnel.tsv",
            qvalues=f"{CALIBRATION_DIR}/qvalues.tsv",
            score_dists=f"{CALIBRATION_DIR}/score_dists.tsv",
            summary=f"{CALIBRATION_DIR}/summary.json"
        log:
            "logs/calibrate.log"
        params:
            thresholds=CALIBRATION.as_dict(),
            two_stage=TWO_STAGE
        threads:
            1
        run:
            os.makedirs(os.path.dirname(log[0]), exist_ok=True)
            arm_inputs = [
                f"{arm}:{wlen}:{path}"
                for (arm, wlen), path in sorted(calibration_arm_inputs().items())
            ]
            cmd = command_tokens("calibration", "python3 -m rnaconsnake.tools.calibration")
            for token in arm_inputs:
                cmd += ["--arm-input", token]
            cmd += [
                "--null-metadata", input.metadata,
                "--output-dir", CALIBRATION_DIR,
                "--rnaz-prob-threshold", str(CALIBRATION.rnaz_prob_threshold),
                "--alifoldz-threshold", str(CALIBRATION.alifoldz_threshold),
                "--rscape-min-pairs", str(CALIBRATION.rscape_min_pairs),
                "--stage1-rnaz-prob", str(CALIBRATION.stage1_rnaz_prob),
                "--locus-min-overlap", str(CALIBRATION.locus_min_overlap),
                "--collapse-ratio-tolerance", str(CALIBRATION.collapse_ratio_tolerance),
                "--dereplicate-method", DEREPLICATE_METHOD,
                "--pair-containment", str(PAIR_CONTAINMENT),
                "--max-container-width", str(CALIBRATION.max_container_width),
                "--container-min-coverage", str(CALIBRATION.container_min_coverage),
                "--representative", CALIBRATION.representative_rule,
            ]
            if TWO_STAGE:
                cmd.append("--two-stage")
            run_checked(cmd, stderr_path=log[0])


rule emit_versions:
    """Record the exact toolchain a run used; a calibrated FDR is only
    reproducible alongside the versions that produced it."""
    output:
        VERSIONS_FILE
    threads:
        1
    run:
        cmd = command_tokens("versions", "python3 -m rnaconsnake.tools.versions") + [
            "--output",
            output[0],
        ]
        for key, command in sorted(SETTINGS.tools.items()):
            cmd += ["--tool", f"{key}={command}"]
        run_checked(cmd)


def benchmark_null_loci():
    """Each null arm's locus tables, one per window length.

    A recovery count on its own is uninterpretable: reported loci cover much of
    an alignment, so the overlap test can be satisfied by chance. The benchmark
    reports how many elements the null arms "recover" too, which needs their
    locus tables here. Grouped per arm, because the tables of one arm pool into
    one baseline sample -- passing them as separate arms would count each window
    length as its own replicate and understate the baseline.
    """
    if not NULL_ARM_ENABLED:
        return {}
    return {
        arm: [
            f"arms/{arm}/generated_files/summary/len_{wlen}/RNAConSnake.nr.csv"
            for wlen in MAXBPSPAN
        ]
        for arm in NULL.null_arms()
    }


rule benchmark_recovery:
    """Positive control: recovery of curated known elements at a given q.

    Not part of ``rule all`` -- it needs a curated truth file and a benchmark
    alignment.  Request it explicitly, or via ``RNAcs --benchmark``.
    """
    input:
        truth=BENCHMARK_TRUTH,
        qvalues=f"{CALIBRATION_DIR}/qvalues.tsv",
        null_loci=[path for paths in benchmark_null_loci().values() for path in paths]
    output:
        f"{BENCHMARK_DIR}/flavivirus_recovery.tsv"
    log:
        "logs/benchmark_recovery.log"
    params:
        alignment=BENCHMARK_ALIGNMENT,
        min_overlap=BENCHMARK_MIN_OVERLAP,
        allow_uncurated=BENCHMARK_ALLOW_UNCURATED,
        null_loci_by_arm=benchmark_null_loci()
    threads:
        1
    run:
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        cmd = command_tokens("benchmark", "python3 -m rnaconsnake.tools.benchmark") + [
            "--truth", input.truth,
            "--qvalues", input.qvalues,
            "--output", output[0],
            "--min-overlap-fraction", str(params.min_overlap),
        ]
        if params.alignment:
            cmd += ["--alignment", str(params.alignment)]
        if params.allow_uncurated:
            cmd.append("--allow-uncurated")
        for arm in sorted(params.null_loci_by_arm):
            cmd += ["--null-loci", *params.null_loci_by_arm[arm]]
        run_checked(cmd, stderr_path=log[0])


rule alignment_report:
    """Which regions of the input alignment can be screened at all.

    Where most sequences are gaps there is no consensus to fold, so no
    candidate can ever be reported -- indistinguishable, in the output, from a
    genuine absence of structure. Stating it explicitly is the difference
    between "nothing is there" and "we could not look".
    """
    input:
        initial_alignment_input
    output:
        report=A("generated_files/alignment_screenability.tsv"),
        metadata=A("generated_files/alignment_screenability.json")
    threads:
        1
    run:
        run_checked(
            command_tokens("alignment_report", "python3 -m rnaconsnake.tools.alignment_report")
            + [
                "--alignment", input[0],
                "--output", output.report,
                "--metadata", output.metadata,
            ]
        )


rule dereplicate_summary:
    """Collapse overlapping windows into one representative per locus.

    RNALalifold reports a real element repeatedly -- near its true extent, and
    again as shorter windows over its individual stable helices, which score
    well in their own right. This keeps the full per-window table untouched and
    adds a non-redundant one alongside it, recording every collapsed member so
    nothing is silently dropped.
    """
    input:
        A("generated_files/summary/len_{wlen}/RNAConSnake.log.csv")
    output:
        nr=A("generated_files/summary/len_{wlen}/RNAConSnake.nr.csv"),
        metadata=A("generated_files/summary/len_{wlen}/RNAConSnake.nr.json")
    log:
        A("generated_files/summary/len_{wlen}/dereplicate.log")
    params:
        # Declared so Snakemake's params rerun-trigger sees a changed method;
        # values read straight from globals inside run: are invisible to it.
        method=DEREPLICATE_METHOD,
        pair_containment=PAIR_CONTAINMENT,
        min_overlap=CALIBRATION.locus_min_overlap,
        max_container_width=CALIBRATION.max_container_width,
        container_min_coverage=CALIBRATION.container_min_coverage,
        representative=CALIBRATION.representative_rule
    threads:
        1
    run:
        run_checked(
            command_tokens("dereplicate", "python3 -m rnaconsnake.tools.dereplicate")
            + [
                "--input", input[0],
                "--output", output.nr,
                "--metadata", output.metadata,
                "--method", params.method,
                "--pair-containment", str(params.pair_containment),
                "--min-overlap", str(params.min_overlap),
                "--max-container-width", str(params.max_container_width),
                "--container-min-coverage", str(params.container_min_coverage),
                "--representative", str(params.representative),
                "--label", f"len{wildcards.wlen}",
            ],
            stdout_path=log[0],
        )


rule render_summary_markdown:
    """Human-readable report: non-redundant loci first, then every window.

    Downstream of de-replication, because the leading block is the collapsed
    locus table. The per-window table follows it unchanged.
    """
    input:
        full=A("generated_files/summary/len_{wlen}/RNAConSnake.log.csv"),
        nr=A("generated_files/summary/len_{wlen}/RNAConSnake.nr.csv")
    output:
        A("generated_files/summary/len_{wlen}/RNAConSnake.md")
    params:
        method=DEREPLICATE_METHOD
    threads:
        1
    run:
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "render-markdown",
                "--label", f"len_{wildcards.wlen}",
                "--nr", input.nr,
                "--full", input.full,
                "--output", output[0],
                "--method", params.method,
            ]
        )


rule clean:
    message: "removing directories: {params}"
    params:
        "Lalifold",
        "generated_files",
        "arms",
        "null_pool",
        "results",
    shell:
        "rm -rf {params}"
