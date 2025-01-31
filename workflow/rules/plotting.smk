rule plot_mapping_rates:
    input:
        expand("results/mapped/{sample}/{sample}Log.final.out", sample=SAMPLES)
    output:
        report("results/plots/mapping_rates.pdf", caption="report/mapping_rates.rst", category="Mapping rates")
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/mapping_rates.log"
    script:
        "../scripts/mapping_rates.R"


rule plot_pca:
    input:
        "results/deseq2/dds.RData",
    output:
        report("results/plots/pca.pdf", caption="report/pca.rst", category="PCA"),
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/pca.log"
    script:
        "../scripts/pca.R"


rule plot_sample_distance:
    input:
        "results/deseq2/dds.RData",
    output:
        report("results/plots/sample_distance.pdf", caption="report/sample_distance.rst", category="Sample distances"),
    params:
        genome=resources.genome,
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/sample_distance.log"
    script:
        "../scripts/heatmap_sd.R"


rule plot_volcano:
    input:
        csv="results/deseq2/{comparison}_{type}.csv",
    output:
        pdf=report("results/plots/volcano/{comparison}_{type}.pdf", caption="report/volcano.rst", category="Volcano plots"),
    params:
        fdr=config["fdr_cutoff"],
        fc=config["fc_cutoff"],
    conda:
        "../envs/deseq2.yml"
    threads: 1
    resources:
        runtime=10
    log:
        "logs/plots/volcano_{comparison}_{type}.log"
    script:
        "../scripts/volcano.R"   


rule plot_te_classes:
    input:
        te_csv=expand("results/deseq2/{comparison}_te.csv", comparison=COMPARISONS)
    output:
        pdf=report(expand("results/plots/te_classes/{comparison}.pdf", comparison=COMPARISONS), caption="report/te_classes.rst", category="TE classes"),
    params:
        fdr=config["fdr_cutoff"],
        lfc=config["fc_cutoff"]
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/te_classes.log"
    script:
        "../scripts/te_classes.R"

