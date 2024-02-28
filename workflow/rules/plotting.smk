rule mapping_rates_plot:
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


rule pca_plot:
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


rule heatmap_sample_distance:
    input:
        "results/deseq2/dds.RData",
    output:
        report("results/plots/sample_distance.pdf", caption="report/sample_distance.rst", category="Sample distances"),
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/sample_distance.log"
    script:
        "../scripts/heatmap_sd.R"


rule volcano_plot:
    input:
        genes="results/deseq2/deseq2_genes.xlsx",
        te="results/deseq2/deseq2_te.xlsx",
    output:
        genes=report(directory("results/plots/volcano_genes/"), caption="report/volcano.rst", category="Volcano plots for genes"),
        te=report(directory("results/plots/volcano_te/"), caption="report/volcano.rst", category="Volcano plots for TEs"),
    params:
        fdr=config["fdr_cutoff"],
        fc=config["fc_cutoff"]
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/volcano.log"
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

