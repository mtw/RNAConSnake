import json
import os
import shlex
import shutil
import subprocess
from pathlib import Path

from snakemake.io import expand


configfile: "config.yaml"

INPUT_ALIGNMENT = config.get("input_alignment")
MAXBPSPAN = config.get("maxbpspan", [100, 200])
LALIFOLD_THREADS = config.get("lalifold_threads", 1)
REMOVE_GAPONLY_GAPRATIO = float(config.get("remove_gaponly_gapratio", 0.5))
REMOVE_GAPONLY_MAX_N = int(config.get("remove_gaponly_max_n", 0))
DO_CM = config.get("do_cm", False)
DO_LOCARNATE = config.get("do_locarnate", False)
DO_PNG = config.get("do_png", True)
DO_RSCAPE = config.get("do_rscape", False)
CM_RNAZ_THRESHOLD = float(config.get("cm_rnaz_prob_threshold", 0.9))
CM_ALIFOLDZ_THRESHOLD = float(config.get("cm_alifoldz_threshold", -2.0))
TOOLS = config.get("tools", {})

wildcard_constraints:
    wlen=r"\d+"


def command_tokens(name, default):
    return shlex.split(str(TOOLS.get(name, default)))


def read_manifest(path):
    with open(path, encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def write_manifest(path, entries):
    with open(path, "w", encoding="utf-8") as handle:
        for entry in entries:
            handle.write(f"{entry}\n")


def write_json(path, payload):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def run_checked(cmd, stdin_path=None, stdout_path=None, stderr_path=None, cwd=None):
    stdin_handle = open(stdin_path, encoding="utf-8") if stdin_path else None
    stdout_handle = open(stdout_path, "w", encoding="utf-8") if stdout_path else None
    stderr_handle = open(stderr_path, "w", encoding="utf-8") if stderr_path else None
    try:
        subprocess.run(
            cmd,
            stdin=stdin_handle,
            stdout=stdout_handle,
            stderr=stderr_handle,
            cwd=cwd,
            text=True,
            check=True,
        )
    finally:
        if stdin_handle:
            stdin_handle.close()
        if stdout_handle:
            stdout_handle.close()
        if stderr_handle:
            stderr_handle.close()


def normalize_rnaalifold_side_output(outdir, canonical_path, suffix):
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


def split_file_basenames(wildcards):
    manifest = checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1]
    return [file[:-4] for file in read_manifest(manifest)]


def initial_alignment_input(wildcards):
    if not INPUT_ALIGNMENT:
        raise ValueError("Missing required config value 'input_alignment'. Use rnaconsnake-run --input-alignment /path/to/input.stk")
    return [str(INPUT_ALIGNMENT)]


def orig_outputs(wildcards):
    return expand(
        "generated_files/orig/len_{wlen}/{file}.orig.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def remgap_outputs(wildcards):
    return expand(
        "generated_files/remgap/len_{wlen}/{file}_remgap.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def strip_outputs(wildcards):
    return expand(
        "generated_files/strip/len_{wlen}/{file}_stripped.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def stk_outputs(wildcards):
    return expand(
        "generated_files/stk/len_{wlen}/{file}.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def cm_status_outputs(wildcards):
    return expand(
        "generated_files/cm/len_{wlen}/{file}.cm.status.json",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def summary_json_outputs(wildcards):
    return expand(
        "generated_files/summary/len_{wlen}/{file}.summary.json",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def rscape_outputs(wildcards):
    return expand(
        "generated_files/rscape/len_{wlen}/{file}.rscape.json",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def png_outputs(wildcards):
    return expand(
        [
            "generated_files/png/len_{wlen}/{file}_aln.png",
            "generated_files/png/len_{wlen}/{file}_ss.png",
        ],
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
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
        expand("generated_files/summary/len_{wlen}/pp_RNALalifold.log", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/pp_RNALalifold.log.csv", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/pp_RNALalifold.md", wlen=MAXBPSPAN),
        expand("generated_files/summary/len_{wlen}/pp_RNALalifold.html", wlen=MAXBPSPAN),
        expand("generated_files/png/len_{wlen}/manifest.txt", wlen=MAXBPSPAN) if DO_PNG else []


rule RNALalifold:
    input:
        initial_alignment_input
    output:
        stdout="Lalifold/len_{wlen}/RNALalifold.out",
        stderr="Lalifold/len_{wlen}/RNALalifold.err",
        multistk="Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    params:
        cmd=TOOLS.get("rnalalifold", "RNALalifold"),
        input_abs=lambda wildcards, input: os.path.abspath(input[0])
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
            -f S \
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
        Path(output.orig).parent.mkdir(parents=True, exist_ok=True)
        Path(output.remgap).parent.mkdir(parents=True, exist_ok=True)
        Path(output.strip).parent.mkdir(parents=True, exist_ok=True)
        Path(output.stk).parent.mkdir(parents=True, exist_ok=True)

        run_checked(["cp", input[0], output.orig])
        run_checked(
            command_tokens("remove_gaponly", "python3 -m rnaconsnake.tools.remove_gaponly")
            + ["-a", input[0], "-i", "stockholm", "-r", str(REMOVE_GAPONLY_GAPRATIO), "-n", str(REMOVE_GAPONLY_MAX_N)],
            stdout_path=output.remgap,
            stderr_path=os.devnull,
        )
        run_checked(
            command_tokens("strip_aln", "python3 -m rnaconsnake.tools.strip_aln")
            + ["-a", output.remgap, "-f", "S", "--nosingle"],
            stdout_path=output.strip,
        )
        run_checked(["cp", output.strip, output.stk])


rule orig_manifest:
    input:
        orig_outputs
    output:
        "generated_files/orig/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


rule remgap_manifest:
    input:
        remgap_outputs
    output:
        "generated_files/remgap/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


rule strip_manifest:
    input:
        strip_outputs
    output:
        "generated_files/strip/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


rule stk_manifest:
    input:
        stk_outputs
    output:
        "generated_files/stk/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


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
        Path(output.aln).parent.mkdir(parents=True, exist_ok=True)
        Path(output.rnaz_txt).parent.mkdir(parents=True, exist_ok=True)
        Path(output.alifoldz_txt).parent.mkdir(parents=True, exist_ok=True)

        run_checked(
            command_tokens("eslreformat", "esl-reformat") + ["clustal", input[0]],
            stdout_path=output.aln,
        )

        cmd = command_tokens("rnaz", "RNAz") + ["-d"]
        if DO_LOCARNATE:
            cmd.append("-l")
        cmd.append(output.aln)
        run_checked(cmd, stdout_path=output.rnaz_txt)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + ["extract-rnaz", "--input", output.rnaz_txt, "--output", output.rnaz_metrics]
        )

        with open(output.aln, encoding="utf-8") as stdin_handle:
            result = subprocess.run(
                command_tokens("alifoldz", "alifoldz.pl") + ["-f", "-t", "0.0"],
                stdin=stdin_handle,
                capture_output=True,
                text=True,
                check=False,
            )

        with open(output.alifoldz_txt, "w", encoding="utf-8") as handle:
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
                + ["extract-alifoldz", "--input", output.alifoldz_txt, "--output", output.alifoldz_metrics]
            )
        else:
            write_json(output.alifoldz_metrics, {"alifoldzscore": ""})


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
        outdir = Path(output.stdout).parent
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
            stdout_path=output.stdout,
            cwd=str(outdir),
        )
        default_stk = outdir / "RNAalifold_results.stk"
        if default_stk.exists() and not Path(output.stk).exists():
            default_stk.rename(output.stk)
        normalize_rnaalifold_side_output(outdir, output.ali_out, "_ali.out")
        normalize_rnaalifold_side_output(outdir, output.dp_ps, "_dp.ps")
        normalize_rnaalifold_side_output(outdir, output.aln_ps, "_aln.ps")
        normalize_rnaalifold_side_output(outdir, output.ss_ps, "_ss.ps")
        run_checked(command_tokens("ps2eps", "ps2eps") + [Path(output.aln_ps).name], cwd=str(outdir))
        run_checked(command_tokens("epstopdf", "epstopdf") + [Path(output.aln_eps).name], cwd=str(outdir))
        run_checked(command_tokens("ps2eps", "ps2eps") + [Path(output.ss_ps).name], cwd=str(outdir))
        run_checked(command_tokens("epstopdf", "epstopdf") + [Path(output.ss_eps).name], cwd=str(outdir))
        for required in [output.stk, output.ali_out, output.dp_ps, output.aln_ps, output.aln_eps, output.aln_pdf, output.ss_ps, output.ss_eps, output.ss_pdf]:
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
        os.makedirs(os.path.dirname(output.aln_png), exist_ok=True)
        run_checked(command_tokens("magick", "magick") + [input.aln_ps, output.aln_png])
        run_checked(command_tokens("magick", "magick") + [input.ss_ps, output.ss_png])


rule png_manifest:
    input:
        png_outputs if DO_PNG else lambda wildcards: []
    output:
        "generated_files/png/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


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
        os.makedirs(os.path.dirname(output.log), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "run-maxcovar",
                "--ali-out",
                input[0],
                "--log",
                output.log,
                "--output",
                output.metrics,
            ]
        )


rule run_rscape_file:
    input:
        "generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        power="generated_files/rscape/len_{wlen}/{file}.power",
        metrics="generated_files/rscape/len_{wlen}/{file}.rscape.json"
    run:
        os.makedirs(os.path.dirname(output.power), exist_ok=True)
        if not DO_RSCAPE:
            Path(output.power).write_text("# R-scape disabled\n", encoding="utf-8")
            write_json(output.metrics, {"rscape_covary_count": ""})
        else:
            outdir = Path(output.power).parent
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
                power_candidates[0].replace(output.power)
                run_checked(
                    command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
                    + [
                        "extract-rscape",
                        "--input",
                        output.power,
                        "--output",
                        output.metrics,
                    ]
                )
            elif len(power_candidates) == 0:
                stdout_text = (workdir / "rscape.stdout").read_text(encoding="utf-8") if (workdir / "rscape.stdout").exists() else ""
                stderr_text = (workdir / "rscape.stderr").read_text(encoding="utf-8") if (workdir / "rscape.stderr").exists() else ""
                Path(output.power).write_text(
                    "# R-scape produced no .power output\n"
                    + f"# exit code: {result.returncode}\n"
                    + (stdout_text if stdout_text else "")
                    + ("\n# stderr\n" + stderr_text if stderr_text else ""),
                    encoding="utf-8",
                )
                write_json(output.metrics, {"rscape_covary_count": "0" if "Number of covarying pairs = 0" in stdout_text else ""})
            else:
                candidates = ", ".join(path.name for path in power_candidates)
                raise FileNotFoundError(f"Could not uniquely identify R-scape .power output in {workdir}: {candidates}")
            sto_pdf = outdir / f"{wildcards.file}.sto.pdf"
            if sto_pdf.exists():
                sto_pdf.unlink()
            metrics = read_json(output.metrics)
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
    input:
        rscape_outputs
    output:
        "generated_files/rscape/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


rule build_cm_file:
    input:
        rnaz="generated_files/rnaz/len_{wlen}/{file}.rnaz.json",
        alifoldz="generated_files/alifoldz/len_{wlen}/{file}.alifoldz.json",
        stk="generated_files/rnaalifold/len_{wlen}/{file}/{file}.RNAalifold_results.stk"
    output:
        "generated_files/cm/len_{wlen}/{file}.cm.status.json"
    run:
        outdir = Path(output[0]).parent
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
                cm_base = str(outdir / wildcards.file)
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

        write_json(output[0], status)


rule cm_manifest:
    input:
        cm_status_outputs
    output:
        "generated_files/cm/len_{wlen}/manifest.txt"
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])


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
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "combine-summary",
                "--wbn",
                wildcards.file,
                "--output",
                output[0],
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
        log="generated_files/summary/len_{wlen}/pp_RNALalifold.log",
        csv="generated_files/summary/len_{wlen}/pp_RNALalifold.log.csv",
        markdown="generated_files/summary/len_{wlen}/pp_RNALalifold.md",
        html="generated_files/summary/len_{wlen}/pp_RNALalifold.html"
    run:
        os.makedirs(os.path.dirname(output.log), exist_ok=True)
        run_checked(
            command_tokens("legacy_postprocess", "python3 -m rnaconsnake.tools.legacy_postprocess")
            + [
                "render-reports",
                "--label",
                f"len_{wildcards.wlen}",
                "--log",
                output.log,
                "--csv",
                output.csv,
                "--markdown",
                output.markdown,
                "--html",
                output.html,
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
