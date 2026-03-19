import os
from snakemake.io import expand

configfile: "config.yaml"
FULLALN = ["aln1"]
MAXBPSPAN = [100,200] #300,400]
LALIFOLD_THREADS = config.get("lalifold_threads", 1)

wildcard_constraints:
    wlen=r"\d+"

rule all:
    input:
#        "DONE",
        expand("Lalifold/len_{wlen}/RC_{wlen}_0001.stk", wlen=MAXBPSPAN),
        expand("Lalifold/len_{wlen}/split/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/remgap/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/stk/len_{wlen}/manifest.txt", wlen=MAXBPSPAN),
        expand("generated_files/strip/len_{wlen}/manifest.txt", wlen=MAXBPSPAN)


#rule filedir:
#    params:
#        lalifold_prefix = config["lalifold_base"],
#        lalifold_stk = "RC_0001.stk",
#        split_dir = "split",

#filedir = rules.filedir.params

rule RNALalifold:
    input:
        expand("data/{sample}.stk",sample=FULLALN)
#    params:
#        EPScols = calc_EPS_cols(config["alen"]),
    output:
        stdout = "Lalifold/len_{wlen}/RNALalifold.out",
        stderr = "Lalifold/len_{wlen}/RNALalifold.err",
        multistk = "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    threads:
        LALIFOLD_THREADS
    shell:
        """
        mkdir -p Lalifold/len_{wildcards.wlen}
        cd Lalifold/len_{wildcards.wlen}
        RNALalifold \
            -L {wildcards.wlen} \
            --aln-stk \
            --id-prefix RC_{wildcards.wlen} \
            --cfactor 0.6 --nfactor 0.5 \
            -r \
            --csv \
            -f S \
            < ../../{input} > RNALalifold.out 2> RNALalifold.err
        """

checkpoint split_stockholm:
    input:
        "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    output:
        directory("Lalifold/len_{wlen}/split"),
        "Lalifold/len_{wlen}/split/manifest.txt"
    threads:
        1
    log:
        out = "Lalifold/len_{wlen}/split/split.out",
        err = "Lalifold/len_{wlen}/split/split.err"
    shell:
        """
        mkdir -p {output[0]}
        cd {output[0]}
        python3 -m rnaconsnake.tools.split_stockholm -a ../RC_{wildcards.wlen}_0001.stk > split.out 2> split.err
        find . -maxdepth 1 -type f -name 'RC_{wildcards.wlen}_*.stk' -print | sed 's#^\./##' | sort > manifest.txt
        """

def read_manifest(path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def write_manifest(path, entries):
    with open(path, "w") as handle:
        for entry in entries:
            handle.write(f"{entry}\n")

def split_file_basenames(wildcards):
    manifest = checkpoints.split_stockholm.get(wlen=wildcards.wlen).output[1]
    return [file[:-4] for file in read_manifest(manifest)]


def remgap_outputs(wildcards):
    return expand(
        "generated_files/remgap/len_{wlen}/{file}_remgap.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def stk_outputs(wildcards):
    return expand(
        "generated_files/stk/len_{wlen}/{file}.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def strip_outputs(wildcards):
    return expand(
        "generated_files/strip/len_{wlen}/{file}_stripped.stk",
        wlen=wildcards.wlen,
        file=split_file_basenames(wildcards),
    )


def call__remove_gaponly(input, output, atype, gapratio):
    command = (
        f"python3 -m rnaconsnake.tools.remove_gaponly -a {input} -i {atype} -r {gapratio} > {output} 2> /dev/null"
    )
    return command

rule remove_gaponly_file:
    input:
        "Lalifold/len_{wlen}/split/{file}.stk"
    output:
        "generated_files/remgap/len_{wlen}/{file}_remgap.stk"
    params:
        atype = "stockholm",
        gapratio = 0.5
    threads:
        1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        shell(call__remove_gaponly(input[0], output[0], params.atype, params.gapratio))

rule remgap_manifest:
    input:
        remgap_outputs
    output:
        "generated_files/remgap/len_{wlen}/manifest.txt"
    threads:
        1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])

rule rename_remgap_file:
    input:
        "generated_files/remgap/len_{wlen}/{file}_remgap.stk"
    output:
        "generated_files/stk/len_{wlen}/{file}.stk"
    threads:
        1
    shell:
        """
        mkdir -p generated_files/stk/len_{wildcards.wlen}
        cp {input} {output}
        """

rule stk_manifest:
    input:
        stk_outputs
    output:
        "generated_files/stk/len_{wlen}/manifest.txt"
    threads:
        1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])

rule strip_alignment_file:
    input:
        "generated_files/remgap/len_{wlen}/{file}_remgap.stk"
    output:
        "generated_files/strip/len_{wlen}/{file}_stripped.stk"
    threads:
        1
    shell:
        """
        mkdir -p generated_files/strip/len_{wildcards.wlen}
        python3 -m rnaconsnake.tools.strip_aln -a {input} -f S --nosingle > {output}
        """

rule strip_manifest:
    input:
        strip_outputs
    output:
        "generated_files/strip/len_{wlen}/manifest.txt"
    threads:
        1
    run:
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        write_manifest(output[0], [os.path.basename(path) for path in input])



# Rule to copy stripped .stk files back to generated_files/stk
#rule copy_strip_to_stk:
#    input:
#    output:
#        expand("generated_files/stk/len_{wlen}", wlen=MAXBPSPAN)
#    run:
####
#            for file in os.listdir(src_dir):
#                if file.endswith(".strip.stk"):
#                    new_file = file.replace(".strip", "")
#                    src = os.path.join(src_dir, file)
##                    shell(f"cp {src} {dest}")

#rule list_files:
#    input:
#        expand("generated_files/remgap/len_{wlen}", wlen=MAXBPSPAN)
#    output:
#        "file_list.txt"
#    run:
#        with open(output[0], 'w') as f:
#            for root, _, files in os.walk("generated_files"):
#                for file in files:
#                    f.write(f"{root}/{file}\n")

#rule create_done:
#    input:
#        "file_list.txt"
#    output:
#        "DONE"
#    shell:
#        """
#        touch {output}
#        """

rule clean:
    message: "removing directories: {params}"
    params:
        "Lalifold",
        "DONE",
        "generated_files",
#        "file_list.txt"
    shell:
        "rm -rf {params}"
