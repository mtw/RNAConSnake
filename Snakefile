configfile: "config.yaml"
FULLALN = ["aln1"]
MAXBPSPAN = [100,200]

#def calc_EPS_cols(alilen):
#    return alilen+10

wildcard_constraints:
    wlen="\d+"

rule all:
#    input:
#        stk=expand("data/{sample}.stk",sample=FULLALN)
    input:
       expand("Lalifold/{wlen}/split/split.out",wlen=MAXBPSPAN),
#        expand("Lalifold/{wlen}/split/RC_0001.stk",wlen=MAXBPSPAN)

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
#    params:
#        EPScols = calc_EPS_cols(config["alen"]),
    output:
        stdout = "Lalifold/{wlen}/RNALalifold.out",
        stderr = "Lalifold/{wlen}/RNALalifold.err",
        multistk = "Lalifold/{wlen}/RC_{wlen}_0001.stk",
    shell:
        """
        RNALalifold                         \
            -L {wildcards.wlen}             \
            --aln-stk                       \
            --id-prefix RC_{wildcards.wlen} \
            --cfactor 0.6 --nfactor 0.5     \
            -r                              \
            --csv                           \
            -f S                            \
            < {input} > {output.stdout} 2> {output.stderr}
         mv RC_{wildcards.wlen}_0001.stk Lalifold/{wildcards.wlen}/
        """

rule move_RNALalifold_data:
    input:
        stdout = "Lalifold/{wlen}/RNALalifold.out",
        stderr = "Lalifold/{wlen}/RNALalifold.err",
    output:
        stdout = "Lalifold/{wlen}/split/RNALalifold.out",
        stderr = "Lalifold/{wlen}/split/RNALalifold.err",
    run:
        shell("mv {input.stdout} {output.stdout}")
        shell("mv {input.stderr} {output.stderr}")

rule move_RNALalifold_stk:
    input:
        "Lalifold/{wlen}/RC_{wlen}_0001.stk"
    output:
        "Lalifold/{wlen}/split/RC_{wlen}_0001.stk"
    run:
       shell("cp {input} {output}")


rule split_stockholm:
    input:
        "Lalifold/{wlen}/split/RC_{wlen}_0001.stk"
#    output:
    #    "Lalifold/{wlen}/split/split.done"
    log:
        out = "Lalifold/{wlen}/split/split.out",
        err = "Lalifold/{wlen}/split/split.err"
    run:
        shell("split_stockholm.pl -a {input} > {log.out} 2> {log.err}")
        shell("mv *stk Lalifold/{wildcards.wlen}/split/")
        shell("rm {input}")
