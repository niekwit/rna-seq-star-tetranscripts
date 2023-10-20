suppressPackageStartupMessages({
  library(pheatmap)
  library(DESeq2)
})

#create output dir
dir.create(file.path(getwd(),"plots-salmon"), showWarnings=FALSE)

print("Generating heatmap count matrix...")

#load DESeq2 data
load(snakemake@input[[1]])

dds <- estimateSizeFactors(dds)
ntd <- normTransform(dds)

#select top 20 genes
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]

#create df for plotting
df <- as.data.frame(colData(dds)[,c("genotype","condition")])

#save heatmap
pdf(snakemake@output[[1]])
pheatmap(assay(ntd)[select,], 
         cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=FALSE, 
         annotation_col=df,
         show_colnames=FALSE)
dev.off()




