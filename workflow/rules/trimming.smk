rule trim_galore:
    input:
        ["reads/{sample}_R1_001.fastq.gz","reads/{sample}_R2_001.fastq.gz"],
    output:
        fasta_fwd=temp("results/trimmed/{sample}_val_1.fq.gz"),
        report_fwd="results/trimmed/{sample}_R1_001.fastq.gz_trimming_report.txt",
        fasta_rev=temp("results/trimmed/{sample}_val_2.fq.gz"),
        report_rev="results/trimmed/{sample}_R2_001.fastq.gz_trimming_report.txt",
    threads: 8
    resources:
        runtime=45,
    log:
        "logs/trim_galore/{sample}.log",
    wrapper:
        "v5.5.1/bio/trim_galore/pe"

        
