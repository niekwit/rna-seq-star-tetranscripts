rule mapping_rates_plot:
    input:
        expand("results/mapped/{sample}/{sample}Log.final.out", sample=SAMPLES)
    output:
        "results/plots/mapping_rates.pdf"
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
        "results/plots/pca.pdf",
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
        "results/plots/sample_distance.pdf",
    conda:
        "../envs/deseq2.yml"
    threads: config["resources"]["plotting"]["cpu"]
    resources:
        runtime=config["resources"]["plotting"]["time"]
    log:
        "logs/plots/sample_distance.log"
    script:
        "../scripts/heatmap_sd.R"

        