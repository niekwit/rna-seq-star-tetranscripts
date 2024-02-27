rule get_readlength:  
    input:
        "results/qc/multiqc/multiqc_data/multiqc_general_stats.txt"
    output:
        "results/qc/readlength.txt",
    conda:
        "../envs/mapping.yml"
    log:
        "logs/mapping/readlength.log"
    script:
        "../scripts/get_readlength.sh"


rule star_index:
    input:
        fasta=resources.fasta,
        gtf=resources.gtf,
        rl="results/qc/readlength.txt",
    output:
        directory(f"resources/{resources.genome}_{resources.build}_index_star/"),
    cache: True
    params:
        sjdbOverhang="$(cat results/qc/readlength.txt)"
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"]
    log:
        "logs/index/star.log"
    wrapper:
        "v3.3.3/bio/star/index"


rule mapping: #change to wrapper later?
    input:
        val1="results/trimmed/{sample}_val_1.fq.gz",
        val2="results/trimmed/{sample}_val_2.fq.gz",
        idx=f"resources/{resources.genome}_{resources.build}_index_star/",
    output:
        "results/mapped/{sample}/{sample}Log.final.out",
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        "results/mapped/{sample}/{sample}ReadsPerGene.out.tab",
    params:
        extra=utils.star_arguments(config)
    threads: config["resources"]["mapping"]["cpu"]
    resources:
        runtime=config["resources"]["mapping"]["time"]
    conda:
        "../envs/mapping.yml"
    log:
        "logs/mapping/{sample}.log"
    shell:
        "rm -rf temp_{wildcards.sample}/ && " # make sure temp dir is not present
        "STAR --runThreadN {threads} "
        "--runMode alignReads "
        "--genomeDir {input.idx} "
        "--readFilesIn {input.val1} {input.val2} "
        "--readFilesCommand zcat "
        "--quantMode TranscriptomeSAM GeneCounts "
        "--twopassMode Basic "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outTmpDir temp_{wildcards.sample}/ "
        "--outFileNamePrefix results/mapped/{wildcards.sample}/{wildcards.sample} "
        "{params.extra} "
        "> {log} 2>&1"


rule index_bam:
    input:
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
    output:
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    threads: config["resources"]["samtools"]["cpu"]
    resources:
        runtime=config["resources"]["samtools"]["time"]
    log:
        "logs/samtools/index_{sample}.log"
    wrapper:
        "v3.3.3/bio/samtools/index"


