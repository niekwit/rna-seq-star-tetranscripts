# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


library(ggrepel)
library(tidyverse)
library(cowplot)

# genes or te
data.type <- snakemake@wildcards[["type"]]

# Load DESeq2 output
csv <- snakemake@input[["csv"]]
df <- read.csv(csv)
contrast <- snakemake@wildcards[["comparison"]]
print(paste("Generating volcano plot for", contrast))

# Get plotting paramaters
fdr <- -log10(snakemake@params[["fdr"]])
lfc <- snakemake@params[["fc"]]

# Log transform padj
df$log.padj <- -log10(df$padj)

# Remove genes with NA in log.padj
df <- df[!(is.na(df$log.padj)), ]

# add fill colour based on log2FC and pvalue
df <- df %>%
  mutate(colour = case_when(
    log2FoldChange > lfc & log.padj > fdr ~ "red",
    log2FoldChange < -lfc & log.padj > fdr ~ "navy",
    log2FoldChange < lfc & log.padj < fdr ~ "grey40",
    log2FoldChange > -lfc & log.padj < fdr ~ "grey40",
  )) 

# select top 8 down and up regulated for labels
df.up <- df %>%
  filter(log2FoldChange > lfc) %>%
  filter(log.padj > fdr) %>%
  arrange(desc(log.padj))
if (nrow(df.up) > 8){
  # Select top 8 rows
  df.up <- df.up[1:8,]
}
  
df.down <- df %>%
  filter(log2FoldChange < -lfc) %>%
  filter(log.padj > fdr) %>%
  arrange(desc(log.padj))
if (nrow(df.down) > 8){
  # Select top 8 rows
  df.down <- df.down[1:8,]
}

df.label <- rbind(df.up, df.down)

# create plot
p <- ggplot(df, aes(x = `log2FoldChange`,
               y = `log.padj`)
       ) + 
  theme_cowplot() +
  theme(axis.text=element_text(size = 16),
        axis.title=element_text(size = 16),
        plot.title = element_text(hjust = 0.5,
                                  size = 16),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        strip.text.x = element_text(size = 16),
        panel.spacing.x = unit(2, "lines")) +
  geom_point(alpha = 0.5,
             shape = 21,
             size = 5,
             colour = "black",
             fill = df$colour) +
  ylab("-log10(adj. p value)") +
  geom_vline(xintercept = c(-lfc, lfc),
             linetype = "dashed", 
             color = "red", 
             linewidth = 0.5) +
  geom_hline(yintercept = fdr,
             linetype = "dashed", 
             color = "red", 
             linewidth = 0.5) +
  ggtitle(paste0(contrast, ": ", data.type))

# Add labels to plot
if (data.type == "te"){
  label_name <- "ensembl_gene_id"
} else {
  label_name <- "external_gene_name"
}
p <- p + geom_label_repel(size = 5,
                   aes(x = `log2FoldChange`,
                       y = `log.padj`,
                       label = .data[[label_name]]), 
                   data = df.label,
                   nudge_x = -0.125,
                   nudge_y = 0.05) +
  scale_fill_manual(values = df$`colour`)
      
# Save plot to file
ggsave(snakemake@output[["pdf"]], p)


sink(log, type = "output")
sink(log, type = "message")
