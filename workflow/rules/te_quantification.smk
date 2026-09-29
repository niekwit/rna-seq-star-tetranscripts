if not config["split_bam"]:

    rule TE_count:
        input:
            bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
            gtf=resources.gtf,
            te_gtf=resources.tegtf,
        params:
            strand=config["strand"],
        output:
            "results/te_count/{sample}.cntTable",
        threads: 3
        resources:
            runtime=165,
        conda:
            "../envs/te.yml"
        log:
            "logs/te_count/{sample}.log",
        shell:
            "TEcount --BAM {input.bam} "
            "--GTF {input.gtf} "
            "--TE {input.te_gtf} "
            "--stranded {params.strand} "
            "--sortByPos "
            "--project {wildcards.sample} "
            "--outdir results/te_count/ "
            "> {log} 2>&1"

else:

    rule split_bam:
        input:
            bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        output:
            bam=temp(
                "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out_chr_{chr}.bam"
            ),
            idx=temp(
                "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out_chr_{chr}.bam.bai"
            ),
        params:
            region=lambda wildcards: wildcards.chr,
            extra="-b",
        threads: 2
        resources:
            runtime=5,
        wrapper:
            "v5.5.1/bio/samtools/view"

    rule TE_count:
        input:
            bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out_chr_{chr}.bam",
            bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out_chr_{chr}.bam.bai",
            gtf=index_resource("gtf"),
            te_gtf=resources.tegtf,
        params:
            strand=config["strand"],
        output:
            temp("results/te_count/{sample}_chr_{chr}.cntTable"),
        threads: 3
        resources:
            runtime=165,
        conda:
            "../envs/te.yml"
        log:
            "logs/te_count/{sample}_chr_{chr}.log",
        shell:
            "TEcount --BAM {input.bam} "
            "--GTF {input.gtf} "
            "--TE {input.te_gtf} "
            "--stranded {params.strand} "
            "--sortByPos "
            "--project {wildcards.sample}_chr_{wildcards.chr} "
            "--outdir results/te_count/ "
            "> {log} 2>&1"

    rule merge_TE_counts:
        input:
            expand("results/te_count/{{sample}}_chr_{chr}.cntTable", chr=CHROMOSOMES),
        output:
            "results/te_count/{sample}.cntTable",
        conda:
            "../envs/te.yml"
        log:
            "logs/te_count/merge_{sample}.log",
        threads: 1
        resources:
            runtime=5,
        shell:
            # Merge all chromosome counts into one file
            # Remove header from all files
            # Sort
            # Add header to first line of concatenated file
            r"cat {input} | "
            "grep -v 'gene/TE' | "
            "sort | "
            "sed '1i gene/TE\tcount' > {output} 2> {log}"


if config["telocal"]["apply"]:

    rule TElocal_count:
        input:
            bam="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
            bai="results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
            gtf=resources.gtf,
            te_index=telocal_index(),
        params:
            strand=config["strand"],
            extra=config["telocal"]["extra_params"],
        output:
            "results/telocal_count/{sample}.cntTable",
        threads: 3
        resources:
            runtime=180,
        conda:
            "../envs/te.yml"
        log:
            "logs/telocal_count/{sample}.log",
        shell:
            # TElocal has no --outdir option, so the output directory is baked into --project
            "TElocal --BAM {input.bam} "
            "--GTF {input.gtf} "
            "--TE {input.te_index} "
            "--stranded {params.strand} "
            "--sortByPos "
            "--project results/telocal_count/{wildcards.sample} "
            "{params.extra} "
            "> {log} 2>&1"
