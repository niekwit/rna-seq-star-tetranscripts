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
        runtime=config["resources"]["fastqc"]["time"],
        mem_mb = 1024
    wrapper:
        "v3.3.3/bio/fastqc"


rule multiqc: 
    input:
        expand("results/qc/fastqc/{sample}{end}_fastqc.zip", sample=SAMPLES, end=["_R1_001","_R2_001"])
    output:
        r="results/qc/multiqc/multiqc.html",
        d=directory("results/qc/multiqc/"),
        t="results/qc/multiqc/multiqc_data/multiqc_general_stats.txt",
    params:
        extra="",  # Optional: extra parameters for multiqc
    log:
        "logs/multiqc/multiqc.log"
    threads: config["resources"]["fastqc"]["cpu"]
    resources:
        runtime=config["resources"]["fastqc"]["time"]
    conda:
        "../envs/mapping.yml"
    shell:
        "multiqc " 
        "--force "
        "--outdir {output.d} "
        "-n multiqc.html "
        "{params.extra} "
        "{input} "
        "> {log} 2>&1"
        