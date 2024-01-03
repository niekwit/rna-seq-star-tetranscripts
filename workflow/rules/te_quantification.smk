rule TE_count:
    input:
        bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
        gtf=resources.gtf,
        te_gtf=resources.tegtf,
    params:
        strand=config["strand"]
    output:
        "results/te_count/{sample}.cntTable",
    threads:
        config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"]
    conda:
        "../envs/te.yml"
    log:
        "logs/te_count/{sample}.log"
    shell:
        "TEcount --BAM {input.bam} "
        "--GTF {input.gtf} "
        "--TE {input.te_gtf} "
        "--stranded {params.strand} "
        "--sortByPos "
        "--project {wildcards.sample} "
        "--outdir results/te_count/ "
        "> {log} 2>&1"

