rule fastqc:
    input:
        "reads/{sample}{end}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}{end}.html",
        zip="results/qc/fastqc/{sample}{end}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "logs/fastqc/{sample}{end}.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"]
    wrapper:
        "v2.0.0/bio/fastqc"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1_001","_R2_001"])
    output:
        "results/qc/multiqc.html",
        "results/qc/multiqc_data/multiqc_general_stats.txt"
    params:
        extra="",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/mapping.yml"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"]
    shell:
        "multiqc " 
        "--force "
        "-o results/qc "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"
