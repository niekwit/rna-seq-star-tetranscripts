rule bam_to_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        idx="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        binsize=config["deeptools"]["binsize"],
        norm=config["deeptools"]["normalisation"],
    threads: config["resources"]["deeptools"]["cpu"]
    resources:
        runtime=config["resources"]["deeptools"]["time"],
    log:
        "logs/deeptools/{sample}.log",
    conda:
        "../envs/deeptools.yml"
    shell:
        "bamCoverage "
        "-p {threads} "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.binsize} "
        "--normalizeUsing {params.norm} "
        "2> {log}"
