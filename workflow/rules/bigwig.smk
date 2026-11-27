rule bam_to_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        binsize=config["deeptools"]["binsize"],
        norm=config["deeptools"]["normalisation"],
        extra=config["deeptools"].get("extra_params", ""),
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
        "{params.extra} "
        "2> {log}"


rule bam_to_stranded_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    output:
        plus="results/bigwig/stranded/{sample}_plus.bw",
        minus="results/bigwig/stranded/{sample}_minus.bw",
    params:
        binsize=config["deeptools"]["binsize"],
        norm=config["deeptools"]["normalisation"],
        extra=config["deeptools"].get("extra_params", ""),
        strand=config["strand"],
    threads: 4
    resources:
        runtime=60,
    log:
        "logs/deeptools/{sample}_stranded.log",
    conda:
        "../envs/deeptools.yml"
    script:
        "../scripts/stranded_bigwig.py"
