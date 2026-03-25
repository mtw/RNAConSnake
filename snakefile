import os
import shutil
import subprocess
from pathlib import Path

from rnaconsnake.workflow_helpers import (
    WorkflowSettings,
    CandidatePaths,
    candidate_outputs_for_manifest,
    initial_alignment_format_code,
    initial_alignment_input as required_initial_alignment_input,
    normalize_rnaalifold_side_output,
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

wildcard_constraints:
    wlen=r"\d+"


def command_tokens(name, default):
    return SETTINGS.command_tokens(name, default)


def write_output_manifest(output_path, input_paths):
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    write_manifest(output_path, [os.path.basename(p) for p in input_paths])


def split_file_basenames(wildcards):
    manifest = checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1]
    return split_file_basenames_from_manifest(manifest)


def initial_alignment_input(wildcards):
    return required_initial_alignment_input(INPUT_ALIGNMENT)


def orig_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.orig,
    )


def remgap_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.remgap,
    )


def strip_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.strip,
    )


def stk_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.stk,
    )


def cm_status_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.cm_status_json,
    )


def summary_json_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.summary_json,
    )


def rscape_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: paths.rscape_json,
    )


def png_outputs(wildcards):
    return candidate_outputs_for_manifest(
        checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1],
        wildcards.wlen,
        lambda paths: [paths.png_aln, paths.png_ss],
    )


rule all:
    input:
        expand("Lalifold/len_{wlen}/RC_{wlen}_0001.stk", wlen=MAXBPSPAN),
        expand("Lalifold/len_{wlen}/split/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/orig/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/remgap/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/strip/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/stk/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/rscape/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/cm/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/RNAConSnake.log", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/RNAConSnake.log.csv", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/RNAConSnake.md", wlen=MAXBPSPAN),
        expand("generated_files/png/len_{wlen}/manifest.txt", wlen=MAXBPSPAN) if DO_PNG else []


rule RNALalifold:
    input:
        initial_alignment_input
    output:
        stdout="Lalifold/len_{wlen}/RNALalifold.out",
        stderr="Lalifold/len_{wlen}/RNALalifold.err",
        multistk="Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    params:
        cmd=SETTINGS.tools.get("rnalalifold", "RNALalifold"),
        input_abs=lambda wildcards, input: os.path.abspath(input[0]),
        input_format=lambda wildcards: INPUT_ALIGNMENT_FORMAT
    threads:
        LALIFOLD_THREADS
    shell:
        """
        mkdir -p Lalifold/len_{wildcards.wlen}
        cd Lalifold/len_{wildcards.wlen}
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
        "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    output:
        directory("Lalifold/len_{wlen}/split"),
        "Lalifold/len_{wlen}/split/manifest.txt"
    log:
        out="Lalifold/len_{wlen}/split/split.out",
        err="Lalifold/len_{wlen}/split/split.err"
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
        "Lalifold/len_{wlen}/split/{file}.stk"
    output:
        orig="generated_files/orig/len_{wlen}/{file}.orig.stk",
        remgap="generated_files/remgap/len_{wlen}/{file}_remgap.stk",
        strip="generated_files/strip/len_{wlen}/{file}_stripped.stk",
        stk="generated_files/stk/len_{wlen}/{file}.stk"
    threads:
        1
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
    output: "generated_files/orig/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule remgap_manifest:
    input: remgap_outputs
    output: "generated_files/remgap/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule strip_manifest:
    input: strip_outputs
    output: "generated_files/strip/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule stk_manifest:
    input: stk_outputs
    output: "generated_files/stk/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule analyze_alignment_file:
    input:
        "generated_files/stk/len_{wlen}/{file}.stk"
    output:
        aln="generated_files/aln/len_{wlen}/{file}.aln",
        rnaz_txt="generated_files/rnaz/len_{wlen}/{file}.rnaz.txt",
        rnaz_metrics="generated_files/rnaz/len_{wlen}/{file}.rnaz.json",
        alifoldz_txt="generated_files/alifoldz/len_{wlen}/{file}.alifoldz.txt",
        alifoldz_metrics="generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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

        with open(paths.aln, encoding="utf-8") as stdin_handle:
            result = subprocess.run(
                command_tokens("alifoldz", "alifoldz.pl") + ["-f", "-t", "0.0"],
                stdin=stdin_handle,
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
            write_json(paths.alifoldz_json, {"alifoldzscore": "0.0"})


rule run_post_rnaalifold_file:
    input:
        "generated_files/stk/len_{wlen}/{file}.stk"
    output:
        stdout="generated_files/rnaalifold/len_{wlen}/{file}/{file}.alifold.out",
        stk="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk",
        ali_out="generated_files/rnaalifold/len_{wlen}/{file}/{file}_ali.out",
        dp_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_dp.ps",
        aln_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.ps",
        aln_eps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.eps",
        aln_pdf="generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.pdf",
        ss_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.ps",
        ss_eps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.eps",
        ss_pdf="generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.pdf"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
        for stray in [Path("alirna.ps"), outdir / "alirna.ps"]:
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
        aln_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_aln.ps",
        ss_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_ss.ps"
    output:
        aln_png="generated_files/png/len_{wlen}/{file}_aln.png",
        ss_png="generated_files/png/len_{wlen}/{file}_ss.png"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
        os.makedirs(os.path.dirname(paths.png_aln), exist_ok=True)
        run_checked(command_tokens("magick", "magick") + [input.aln_ps, paths.png_aln])
        run_checked(command_tokens("magick", "magick") + [input.ss_ps, paths.png_ss])


rule png_manifest:
    input: png_outputs if DO_PNG else lambda wildcards: []
    output: "generated_files/png/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule reformat_rnaalifold_results_file:
    input:
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln"
    run:
        run_checked(
            command_tokens("eslreformat", "esl-reformat") + ["clustal", input[0]],
            stdout_path=output[0],
        )


rule clean_rnaalifold_clustal_file:
    input:
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln"
    output:
        backup="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.aln~",
        cleaned="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.cleaned.aln"
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
    input:
        aln="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.cleaned.aln",
        dp_ps="generated_files/rnaalifold/len_{wlen}/{file}/{file}_dp.ps"
    output:
        "generated_files/refold/len_{wlen}/{file}_refold.out"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        refold_cmd = command_tokens("refold", "refold.pl") + [input.aln, input.dp_ps]
        rnafold_cmd = command_tokens("rnafold", "RNAfold") + ["--noPS", "-C"]
        with open(output[0], "w", encoding="utf-8") as out_handle:
            first = subprocess.Popen(refold_cmd, stdout=subprocess.PIPE, text=True)
            try:
                subprocess.run(rnafold_cmd, stdin=first.stdout, stdout=out_handle, text=True, check=True)
            finally:
                if first.stdout:
                    first.stdout.close()
                first.wait()
                if first.returncode != 0:
                    raise subprocess.CalledProcessError(first.returncode, refold_cmd)


rule extract_refold_metrics_file:
    input:
        refold="generated_files/refold/len_{wlen}/{file}_refold.out",
        stk="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        "generated_files/refold/len_{wlen}/{file}.refold.json"
    run:
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "extract-refold",
                "--refold-output",
                input.refold,
                "--rnaalifold-stk",
                input.stk,
                "--output",
                output[0],
            ]
        )


rule run_maxcovar_file:
    input:
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}_ali.out"
    output:
        log="generated_files/maxcovar/len_{wlen}/{file}_alifoldmaxcovar.log",
        metrics="generated_files/maxcovar/len_{wlen}/{file}.maxcovar.json"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        power="generated_files/rscape/len_{wlen}/{file}.power",
        metrics="generated_files/rscape/len_{wlen}/{file}.rscape.json"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
    output: "generated_files/rscape/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule build_cm_file:
    input:
        rnaz="generated_files/rnaz/len_{wlen}/{file}.rnaz.json",
        alifoldz="generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json",
        stk="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        "generated_files/cm/len_{wlen}/{file}.cm.status.json"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
    output: "generated_files/cm/len_{wlen}/manifest.txt"
    run: write_output_manifest(output[0], input)


rule combine_summary_metrics_file:
    input:
        rnaz="generated_files/rnaz/len_{wlen}/{file}.rnaz.json",
        alifoldz="generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json",
        refold="generated_files/refold/len_{wlen}/{file}.refold.json",
        maxcov="generated_files/maxcovar/len_{wlen}/{file}.maxcovar.json",
        rscape="generated_files/rscape/len_{wlen}/{file}.rscape.json"
    output:
        "generated_files/summary/len_{wlen}/{file}.summary.json"
    run:
        paths = CandidatePaths(wlen=wildcards.wlen, file=wildcards.file)
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
        log="generated_files/summary/len_{wlen}/RNAConSnake.log",
        csv="generated_files/summary/len_{wlen}/RNAConSnake.log.csv",
        markdown="generated_files/summary/len_{wlen}/RNAConSnake.md"
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
                "--markdown",
                output.markdown,
                *sorted(input),
            ]
        )


rule clean:
    message: "removing directories: {params}"
    params:
        "Lalifold",
        "generated_files",
    shell:
        "rm -rf {params}"
