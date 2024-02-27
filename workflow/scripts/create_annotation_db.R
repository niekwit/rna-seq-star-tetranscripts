# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(rtracklayer)

gtf <- snakemake@input[["gtf"]]

# Gene annotation info
db <- rtracklayer::import(gtf)
gene.info <- data.frame(ensembl_gene_id = db$gene_id, 
                  external_gene_name = db$gene_name, 
                  gene_biotype = db$gene_biotype)

# Save gene annotation info
save(gene.info, file = snakemake@output[["edb"]])

sink(log, type = "output")
sink(log, type = "message")