rule deseq2:
    input:
        counts=expand("results/te_count/{sample}.cntTable", sample=SAMPLES),
        gtf=resources.gtf,
    output:
        genes=report("results/deseq2/deseq2_genes.xlsx", caption="report/deseq2.rst", category="Differential Expression Analysis of genes"),
        te=report("results/deseq2/deseq2_te.xlsx", caption="report/deseq2.rst", category="Differential Expression Analysis of TEs"),
        rdata="results/deseq2/dds.RData",
    retries: 5 # gene annotation may fail due to database failed connection
    params:
        strand=config["strand"],
        genome=resources.genome
    threads: config["resources"]["deseq2"]["cpu"]
    resources:
        runtime=config["resources"]["deseq2"]["time"]
    conda:
        "../envs/deseq2.yml"
    log:
        "logs/deseq2/deseq2.log"
    script:
        "../scripts/deseq2.R"

