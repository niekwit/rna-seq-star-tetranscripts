# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

load(snakemake@input[[1]])

batches <- dds$batch

if (length(unique(batches)) > 1) {
  print("Removing batch effect from data...")
  # Extracting transformed values
  vsd <- vst(dds, blind = FALSE)
  
  # Remove batch variation with limma
  mat <- assay(vsd)
  mm <- model.matrix(~comb, colData(vsd))
  mat <- limma::removeBatchEffect(mat, 
                                  batch = vsd$batch, 
                                  design = mm)
  assay(vsd) <- mat
  
  # Select appropriate colour palette
  if (length(unique(vsd$treatment)) == 2) {
    palette <- "Paired"
  } else {
    palette <- "Dark2"
  }
  
  # Create PCA plot
  p <- plotPCA(vsd, intgroup = c("genotype", "treatment")) +
    geom_text_repel(aes(label = vsd$sample),
                    size = 6) + 
    guides(colour = "none") +
    theme_cowplot(18) +
    scale_color_brewer(palette = palette)

} else {
  print("No correction for batch effect...")
  # Log transform data
  rld <- rlog(dds)

  # Select appropriate colour palette
  conditions <- length(unique(rld$genotype)) * length(unique(rld$treatment))
  if (conditions <= 12) {
    palette <- "Paired"
  } else {
    shapes <- c(16, 15, 17, 18, 1, 0, 2, 5, 6, 3, 4, 7, 8, 9, 10)
    shapes <- shapes[1:conditions]
    palette <- "Dark2"
  }

  # Prepare data for PCA
  # https://github.com/thelovelab/DESeq2/blob/30cc350d58c307a2c564c2b372ec2a50462a8097/R/plots.R#L239
  
  # Calculate the variance for each gene
  rv <- rowVars(assay(rld))
  
  # Select the ntop genes by variance
  ntop <- 500
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(rld)[select,]))
  
  # The contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  # Prepare data frame for plotting
  df <- pca$x[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    left_join(samples, by = "sample")

  p <- ggplot(df, 
              aes(x = PC1, 
                  y = PC2, 
                  color = treatment, 
                  shape = genotype)) +
    geom_point(size = 8) +
    geom_text_repel(aes(label = sample), 
                    size = 4) +
    theme_cowplot(18) +
    scale_color_brewer(palette = palette) +
    scale_shape_manual(values = shapes) +
    labs(x = paste0("PC1 (", round(percentVar[1], 3) * 100, "% variance)"),
         y = paste0("PC2 (", round(percentVar[2], 3) * 100, "% variance)"))
}  

# Save plot to file
ggsave(snakemake@output[[1]], 
       p, 
       width=10,
       height=10)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")