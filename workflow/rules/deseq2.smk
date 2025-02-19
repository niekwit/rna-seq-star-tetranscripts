rule create_annotation_db:
    input:
        gtf=resources.gtf,
    output:
        edb="resources/edb.RData"
    params:
        genome=resources.genome
    conda:
        "../envs/rtracklayer.yml"
    threads: 2
    resources:
        runtime=30
    log:
        "logs/rtracklayer/create_annotation_db.log"
    script:
        "../scripts/create_annotation_db.R"

rule deseq2:
    input:
        counts=expand("results/te_count/{sample}.cntTable", sample=SAMPLES),
        edb="resources/edb.RData",
    output:
        rdata="results/deseq2/dds.RData",
        genes_csv=expand("results/deseq2/{comparison}_genes.csv", comparison=COMPARISONS),
        te_csv=expand("results/deseq2/{comparison}_te.csv", comparison=COMPARISONS)
    params:
        strand=config["strand"],
        genome=resources.genome
    threads: 4
    resources:
        runtime=90
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "../scripts/deseq2.R"

