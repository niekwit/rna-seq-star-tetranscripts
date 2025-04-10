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

if config["spike_in"]["apply"]:
    logger.info("Spike-in applied")
    check_spike_in_resources()
    rule combine_fasta:
        input:
            genome=resources.fasta,
            spikein=config["spike_in"]["fasta"],
        output:
            "resources/combined.fasta",
        log:
            "logs/resources/combine_fasta.log"
        threads: 1
        resources: 
            runtime=5
        conda:
            "../envs/mapping.yml"
        shell:
            "cat {input.genome} {input.spikein} > {output}"


    use rule combine_fasta as combine_gtf with:
            input:
                genome=resources.gtf,
                spikein=config["spike_in"]["gtf"],
            output:
                "resources/combined.gtf",
            log:
                "logs/resources/combine_gtf.log"
