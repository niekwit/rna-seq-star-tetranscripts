rule get_readlength:  
    input:
        r="results/qc/multiqc/multiqc.html",
        d="results/qc/multiqc/",
        t="results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",
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
        fasta=index_resource("fasta"),
        gtf=index_resource("gtf"),
        rl="results/qc/readlength.txt",
    output:
        directory=directory(f"resources/{resources.genome}_{resources.build}_index_star/"),
        files=multiext(f"resources/{resources.genome}_{resources.build}_index_star/", "chrLength.txt", 
        "chrNameLength.txt", 
        "chrName.txt", 
        "chrStart.txt", 
        "exonGeTrInfo.tab", 
        "exonInfo.tab",
        "geneInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"),
    #cache: True
    params:
        sjdbOverhang="$(cat results/qc/readlength.txt)",
        extra="",
    threads: 36
    resources:
        runtime=120,
    conda:
        "../envs/mapping.yml"
    log:
        "logs/index/star.log"
    shell:
        "STAR "
        "--runThreadN {threads} "  # Number of threads
        "--runMode genomeGenerate "  # Indexation mode
        "--genomeFastaFiles {input.fasta} "  # Path to fasta files
        "--sjdbOverhang {params.sjdbOverhang} "  # Read-len - 1
        "--sjdbGTFfile {input.gtf} "
        "{params.extra} "  # Optional parameters
        "--genomeDir {output.directory} "  # Path to output
        "> {log} 2>&1"  # Logging


rule mapping: 
    input:
        val1="results/trimmed/{sample}_val_1.fq.gz",
        val2="results/trimmed/{sample}_val_2.fq.gz",
        idx=f"resources/{resources.genome}_{resources.build}_index_star/",
        idxfiles=multiext(f"resources/{resources.genome}_{resources.build}_index_star/", "chrLength.txt", 
        "chrNameLength.txt", 
        "chrName.txt", 
        "chrStart.txt", 
        "exonGeTrInfo.tab", 
        "exonInfo.tab",
        "geneInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"),
    output:
        "results/mapped/{sample}/{sample}Log.final.out",
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
        "results/mapped/{sample}/{sample}ReadsPerGene.out.tab",
    params:
        prefix=lambda wc, output: output[0].replace("Log.final.out", ""),
        extra=star_arguments(config)
    threads: 12
    resources:
        runtime=90
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
        "--quantMode GeneCounts "
        "--twopassMode Basic "
        "--outSAMunmapped None "
        "--outSAMtype BAM SortedByCoordinate "
        "--outTmpDir temp_{wildcards.sample}/ "
        "--outFileNamePrefix {params.prefix} "
        "{params.extra} "
        "> {log} 2>&1"


rule index_bam:
    input:
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam",
    output:
        "results/mapped/{sample}/{sample}Aligned.sortedByCoord.out.bam.bai",
    threads: 4
    resources:
        runtime=10
    log:
        "logs/samtools/index_{sample}.log"
    wrapper:
        "v5.5.1/bio/samtools/index"

