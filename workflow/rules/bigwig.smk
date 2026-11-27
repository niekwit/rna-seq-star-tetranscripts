rule bam_to_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    output:
        bw="results/bigwig/{sample}.bw",
        # Scale factor computed by bamCoverage, reused for the stranded bigWigs
        scale="results/bigwig/scale_factors/{sample}.txt",
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
        "-o {output.bw} "
        "--binSize {params.binsize} "
        "--normalizeUsing {params.norm} "
        "{params.extra} "
        "--verbose "
        "> {log} 2>&1; "
        "grep -m1 '^Scale factor:' {log} | cut -d' ' -f3 > {output.scale}"


rule bam_to_stranded_bigwig:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
        scale="results/bigwig/scale_factors/{sample}.txt",
    output:
        plus="results/bigwig/stranded/{sample}_plus.bw",
        minus="results/bigwig/stranded/{sample}_minus.bw",
    params:
        binsize=config["deeptools"]["binsize"],
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
