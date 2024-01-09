rule trim_galore:
    input:
        ["reads/{sample}_R1_001.fastq.gz","reads/{sample}_R2_001.fastq.gz"],
    output:
        fasta_fwd=temp("results/trimmed/{sample}_val_1.fq.gz"),
        report_fwd="results/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
        fasta_rev=temp("results/trimmed/{sample}_val_2.fq.gz"),
        report_rev="results/trimmed/{sample}_R2_001.fastq.gz_trimming_report.txt",
    threads: config["resources"]["trim"]["cpu"]
    resources:
        runtime=config["resources"]["trim"]["time"],
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v3.3.3/bio/trim_galore/pe"

        
