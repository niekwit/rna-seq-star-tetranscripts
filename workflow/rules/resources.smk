rule get_fasta:
    output:
        ensure(resources.fasta, sha256=resources.fasta_sha256)
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz 2> {log}"


rule get_gtf: 
    output:
        ensure(resources.gtf, sha256=resources.gtf_sha256)
    retries: 3
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log"
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz 2> {log}"


rule get_te_gtf:
    output:
        ensure(resources.tegtf, sha256=resources.tegtf_sha256)
    retries: 3
    params:
        url=resources.tegtf_url,
    log:
        "logs/resources/get_te_gtf.log"
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz && gunzip -f {output}.gz 2> {log}"


rule combined_gtf:
    input:
        gtf=resources.gtf,
        te_gtf=resources.tegtf,
    output:
        temp("resources/combined.gtf"),
    conda:
        "../envs/mapping.yml"
    log:
        "logs/resources/combined_gtf.log"
    shell:
        "cat {input.gtf} {input.te_gtf} > {output} 2> {log}}"


# much quicker than reading gene names from GTF file in R
rule get_gene_names_from_gtf:
    input:
        "resources/combined.gtf"
    output:
        "resources/all_gene_names.txt"
    threads: 1
    resources:
        runtime=10
    conda:
        "../envs/mapping.yml"
    log:
        "logs/resources/get_gene_names_from_gtf.log"
    shell:
        """
        grep -o -P '(?<=gene_id ")[^"]*' {input} | sort | uniq > {output} 2> {log}
        """

