rule deseq2:
    input:
        counts=expand("results/te_count/{sample}.cntTable", sample=SAMPLES),
        gtf="resources/combined.gtf",
        genes="resources/all_gene_names.txt",
    output:
        genes="results/deseq2/deseq2_genes.xlsx",
        te="results/deseq2/deseq2_te.xlsx",
        rdata="results/deseq2/dds.RData"
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
        "scripts/deseq2.R"

