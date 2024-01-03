rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    threads: 1
    resources: 
        runtime=30
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz 2> {log};"
        "gunzip -f {output}.gz 2>> {log}"


rule get_gtf: 
    output:
        resources.gtf,
    retries: 3
    params:
        url=resources.gtf_url,
    log:
        "logs/resources/get_gtf.log"
    threads: 1
    resources: 
        runtime=30
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz 2> {log};"
        "gunzip -f {output}.gz 2>> {log}"


rule get_te_gtf:
    output:
        resources.tegtf,
    retries: 3
    params:
        url=resources.tegtf_url,
    log:
        "logs/resources/get_te_gtf.log"
    threads: 1
    resources: 
        runtime=30
    conda:
        "../envs/mapping.yml"
    shell:
        "wget -q {params.url} -O {output}.gz 2> {log};"
        "gunzip -f {output}.gz 2>> {log}"


rule compress_resources:
    input:
        f=resources.fasta,
        g=resources.gtf,
        tg=resources.tegtf,
        p="results/plots/pca.pdf",# dummy input to make sure this rule is executed at the end
    output:
        f"{resources.fasta}.gz",
        f"{resources.gtf}.gz",
        f"{resources.tegtf}.gz",
    params:
        pigz_options="",
    threads: config["resources"]["mapping"]["cpu"]
    resources: 
        runtime=config["resources"]["mapping"]["time"]
    log:
        "logs/resources/compress_resources.log"
    conda:
        "../envs/mapping.yml"
    shell:
        "pigz "
        "-p {threads} "
        "{params.pigz_options} "
        "{input.f} "
        "{input.g} "
        "{input.tg} "
        "> {log} 2>&1"