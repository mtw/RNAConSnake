configfile: "config.yaml"
FULLALN = ["aln1"]
MAXBPSPAN = [100]

#def calc_EPS_cols(alilen):
#    return alilen+10

wildcard_constraints:
    wlen="\d+"

rule all:
#    input:
#        stk=expand("data/{sample}.stk",sample=FULLALN)
    input:
        expand("Lalifold/{wlen}/split/RNALalifold.out",wlen=MAXBPSPAN)
        #expand("Lalifold/{wlen}/split/RC_0001.stk",wlen=MAXBPSPAN)

rule filedir:
    params:
        lalifold_prefix = config["lalifold_base"],
        lalifold_stk = "RC_0001.stk",
        split_dir = "split",


filedir = rules.filedir.params

rule RNALalifold:
    input:
        #rules.all.input.stk
        expand("data/{sample}.stk",sample=FULLALN)
    #output:
    #
#    params:
#        maxBPspan="{wildcards.wlen}"
#        EPScols = calc_EPS_cols(config["alen"]),
    output:
        stdout = "Lalifold/{wlen}/RNALalifold.out",
        stderr = "Lalifold/{wlen}/RNALalifold.err",
#multistk = "Lalifold/{wlen}/RC_0001.stk",
    run:
        shell("RNALalifold -L {wildcards.wlen} "
        "--aln-stk "
        "--id-prefix=RC "
        "--cfactor 0.6 --nfactor 0.5 "
        "-r "
        "--csv "
        "-f S "
        "< {input} > {output.stdout} 2> {output.stderr}")
    #    shell("mv RC_0001.stk Lalifold/{wlen}/")


rule move_RNALalifold_data:
    input:
        stdout = "Lalifold/{wlen}/RNALalifold.out",
        stderr = "Lalifold/{wlen}/RNALalifold.err",
    output:
        stdout = "Lalifold/{wlen}/split/RNALalifold.out",
        stderr = "Lalifold/{wlen}/split/RNALalifold.err",
    run:
        shell("mv {input.stdout} Lalifold/{wlen}/split/")
        shell("mv {input.stderr} Lalifold/{wlen}/split/")

rule move_RNALalifold_stk:
    input:
        multistk = "RC_0001.stk",
    output:
        multistk = "RC_0001.stk",
    run:
       shell("mv {input.multistk} Lalifold/{wlen}/split/")


#rule split_stockholm:
#    input:
#        rules.move_RNALalifold_data.output
#    output:
#        "{filedir.lalifold_prefix}/{wlen}/{filedir.split_dir}/split.done"
#    log:
#        out = "split.out",
#        err = "split.err"
#    run:
    #    shell("split_stockholm.pl -a {input} > {log.out} 2> {log.err}")
    #    shell("mv *stk Lalifold/{wlen}/split/")
#        shell("touch {output}")
