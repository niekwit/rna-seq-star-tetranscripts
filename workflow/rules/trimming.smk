rule trim_galore:
    input:
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz",
    output:
        temp("results/trimmed/{sample}_val_1.fq.gz"),
        "results/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
        temp("results/trimmed/{sample}_val_2.fq.gz"),
        "results/trimmed/{sample}_R2_001.fastq.gz_trimming_report.txt",
    threads: config["resources"]["trim"]["cpu"]
    conda:
        "../envs/mapping.yml"
    log:
        "logs/trim_galore/{sample}.log",
    shell:
        "trim_galore -j {threads} "
        "-q 20 " 
        "--basename {wildcards.sample} "
        "-o results/trimmed "
        "--paired {input.r1} {input.r2} " 
        "> {log} 2>&1"

        
