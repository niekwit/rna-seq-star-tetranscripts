# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

# dds variable is loaded from deseq2.R
load(snakemake@input[[1]])

# Get batch information
samples <- read.csv("config/samples.csv", header = TRUE)
batches <- unique(samples$batch)

# Extracting transformed values
vsd <- vst(dds, blind = FALSE)

# Remove batch effect if multiple batches are present
if (length(unique(batches)) > 1) {
  print("Removing batch effect from data...")
  mat <- assay(vsd)
  mm <- model.matrix(~comb, colData(vsd))
  mat <- limma::removeBatchEffect(mat, 
                                  batch = vsd$batch, 
                                  design = mm)
  assay(vsd) <- mat
} else {print("No correction for batch effect...")}

# Select colours and shapes for plotting
shapes <- c(16, 15, 17, 18, 1, 0, 2, 5, 6, 3, 4, 7, 8, 9, 10)
shapes <- shapes[1:length(unique(vsd$genotype))]
if (length(unique(vsd$treatment)) <= 8) {
  palette <- "Dark2"
} else {
  palette <- "Paired"
}

# Prepare data for PCA
# https://github.com/thelovelab/DESeq2/blob/30cc350d58c307a2c564c2b372ec2a50462a8097/R/plots.R#L239

# Calculate the variance for each gene
rv <- rowVars(assay(vsd))

# Select the ntop genes by variance
ntop <- 500
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

# Perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))

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

# Save plot to file
ggsave(snakemake@output[[1]], 
       p, 
       width=10,
       height=10)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")