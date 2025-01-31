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
        runtime=15
    conda:
        "../envs/mapping.yml"
    script:
        "../scripts/get_resource.sh"


use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


use rule get_fasta as get_te_gtf with:
        output:
            resources.tegtf,
        params:
            url=resources.tegtf_url,
        log:
            "logs/resources/get_te_gtf.log"


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
    threads: 4
    resources: 
        runtime=15
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
        