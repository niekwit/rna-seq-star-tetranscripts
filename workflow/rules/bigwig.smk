rule bam_to_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        binsize=config["deeptools"]["binsize"],
        norm=config["deeptools"]["normalisation"],
    threads: 4
    resources:
        runtime=30,
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
