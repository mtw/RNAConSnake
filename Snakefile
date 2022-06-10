configfile: "config.yaml"
FULLALN = ["aln1"]
MAXBPSPAN = [100,150]

#def calc_EPS_cols(alilen):
#    return alilen+10

rule all:
    input:
        expand("Lalifold/{wlen}/split/split.done",wlen=MAXBPSPAN)

rule filedir:
    params:
        lalifold_prefix = config["lalifold_base"],
        split_dir = "split"

filedir = rules.filedir.params

rule RNALalifold:
    input:
        stk=expand("data/{sample}.stk",sample=FULLALN)
    output:
        multistk="RC_0001.stk"
    params:
        alen = {wlen}
#        alen = config["alen"],
#        EPScols = calc_EPS_cols(config["alen"]),
    log:
        out = "RNALalifold.out",
        err = "RNALalifold.err"
    shell:
        "RNALalifold -L {params.alen} "
        "--aln-stk "
        "--id-prefix=RC "
        "--cfactor 0.6 --nfactor 0.5 "
        "-r "
        "--csv "
        "-f S "
        "< {input.stk} > {log.out} 2> {log.err}"

#rule move_RNALalifold_data:
#    input:
#        rules.RNALalifold.output.multistk
#    output:
#        "Lalifold/{wlen}/RC_0001.stk"
#    run:
#        shell("mv {input} {filedir.lalifold_prefix}/{wlen}/")
#        shell("mv {rules.RNALalifold.log.out} {filedir.lalifold_prefix}/{wlen}/")
#        shell("mv {rules.RNALalifold.log.err} {filedir.lalifold_prefix}/{wlen}/")

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
