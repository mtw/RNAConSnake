import os
from snakemake.io import expand, glob_wildcards

configfile: "config.yaml"
FULLALN = ["aln1"]
MAXBPSPAN = [100,200] #300,400]

wildcard_constraints:
    wlen=r"\d+"

rule all:
    input:
#        "DONE",
        expand("Lalifold/len_{wlen}/RC_{wlen}_0001.stk", wlen=MAXBPSPAN),
        expand("Lalifold/len_{wlen}/split/split.done", wlen=MAXBPSPAN),
        expand("generated_files/stk/len_{wlen}", wlen=MAXBPSPAN),
        expand("generated_files/strip/len_{wlen}", wlen=MAXBPSPAN)


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
        #multistk = "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
        multistk = "RC_{wlen}_0001.stk"
    shell:
        """
        RNALalifold \
            -L {wildcards.wlen} \
            --aln-stk \
            --id-prefix RC_{wildcards.wlen} \
            --cfactor 0.6 --nfactor 0.5 \
            -r \
            --csv \
            -f S \
            < {input} > {output.stdout} 2> {output.stderr}
        """

rule movestk:
    input:
        rules.RNALalifold.output.multistk
    output:
        "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    shell:
       "mv {input} {output}"

rule split_stockholm:
    input:
        "Lalifold/len_{wlen}/RC_{wlen}_0001.stk"
    output:
        directory("Lalifold/len_{wlen}/split"),
        touch("Lalifold/len_{wlen}/split/split.done")
    log:
        out = "Lalifold/len_{wlen}/split/split.out",
        err = "Lalifold/len_{wlen}/split/split.err"
    shell:
        """
        mkdir -p {output[0]}
        split_stockholm.pl -a {input} > {log.out} 2> {log.err}
        mv RC_{wildcards.wlen}_*.stk {output[0]}/
        """
#    run:
#        shell("split_stockholm.pl -a {input} > {log.out} 2> {log.err}")
#        shell("touch split.done")
#        shell("mv split.done Lalifold/len_{wildcards.wlen}/split/")
#        #shell("rm {input}")
#        shell("mv *stk Lalifold/len_{wildcards.wlen}/split/")

# Use glob_wildcards to find all stk files and group them by wlen
found_files = glob_wildcards("Lalifold/len_{wlen}/split/{file}.stk")

# Create a dictionary to map wlen to their respective files
wlen_to_files = {}
for wlen, file in zip(found_files.wlen, found_files.file):
    wlen = int(wlen)
    if wlen not in wlen_to_files:
        wlen_to_files[wlen] = []
    wlen_to_files[wlen].append(file)

def call__remove_gaponly(input,output,atype,gapratio):
    command = (
        f"remove-gaponly.pl -a {input} -i {atype} -r {gapratio} > {output} 2> /dev/null"
    )
    return command

rule remove_gaponly:
    input:
        "Lalifold/len_{wlen}/split/split.done"
    output:
        directory("generated_files/remgap/len_{wlen}")
    params:
        atype = "stockholm",
        gapratio = 0.5
    run:
        wlen = int(wildcards.wlen)
        output_dir = os.path.join(output[0])
        os.makedirs(output_dir, exist_ok=True)
        for file in wlen_to_files[int(wlen)]:
            inp = f"Lalifold/len_{wlen}/split/{file}.stk"
            out = os.path.join(output_dir, f"{file}.remgap.stk")
            cmd = call__remove_gaponly(inp,out,params.atype,params.gapratio)
            shell(cmd)

rule rename_remgap:
    input:
        expand("generated_files/remgap/len_{wlen}", wlen=MAXBPSPAN)
    output:
        directory("generated_files/stk/len_{wlen}")
    run:
        wlen = int(wildcards.wlen)
        src_dir = f"generated_files/remgap/len_{wlen}"
        dest_dir = f"generated_files/stk/len_{wlen}"
        os.makedirs(dest_dir, exist_ok=True)

        for file in os.listdir(src_dir):
            if file.endswith(".remgap.stk"):
                new_file = file.replace(".remgap", "")
                src = os.path.join(src_dir, file)
                dest = os.path.join(dest_dir, new_file)
                shell(f"cp {src} {dest}")

rule strip_alignment:
    input:
        expand("generated_files/remgap/len_{wlen}", wlen=MAXBPSPAN)
    output:
        directory("generated_files/strip/len_{wlen}")
    run:
        wlen = int(wildcards.wlen)
        input_dir = f"generated_files/remgap/len_{wlen}"
        output_dir = f"generated_files/strip/len_{wlen}"
        stk_dir = f"generated_files/stk/len_{wlen}"
        os.makedirs(output_dir, exist_ok=True)

        for file in os.listdir(input_dir):
            if file.endswith(".stk"):
                new_file = file.replace("remgap", "strip")
                plain_file = file.replace(".remgap.stk", ".stk")
                inp = os.path.join(input_dir, file)
                out = os.path.join(output_dir, new_file)
                stk = os.path.join(stk_dir, plain_file)
                shell(f"strip_aln.pl -a {inp} -f S --nosingle > {out}")
                shell(f"cp -f {out} {stk}")



#rule strip_aln:
#    input:
#        "generated_files/stk/len_{wlen}/{file}.stk"
#    output:
#        "generated_files/strip/len_{wlen}/{file}.strip.stk",
#        directory("generated_files/strip/len_{wlen}")
#    shell:
#        """
#        mkdir -p {output}
#        strip_aln.pl -a {input} > {output}/{wildcards.file}.strip.stk
#        """

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
