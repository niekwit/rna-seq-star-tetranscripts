# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(rtracklayer)
library(dplyr)

gtf <- snakemake@input[["gtf"]]
genome <- snakemake@params[["genome"]]

# Gene annotation info
db <- rtracklayer::import(gtf)

if (genome == "T2T-CHM13v2.0") {
  gene.info <- data.frame(gene_id = db$db_xref,
                          type = db$type,
                          external_gene_name = db$gene, 
                          gene_biotype = db$gene_biotype) %>%
    filter(type == "gene") 
} else {
  gene.info <- data.frame(ensembl_gene_id = db$gene_id,
                          external_gene_name = db$gene_name,
                          gene_biotype = db$gene_biotype,
                          gene_id = db$gene_id)
}



# Save gene annotation info
save(gene.info, file = snakemake@output[["edb"]])

sink(log, type = "output")
sink(log, type = "message")