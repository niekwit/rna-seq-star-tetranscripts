# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(readxl)
library(ggrepel)
library(dplyr)
library(stringr)
library(cowplot)

# get xlsx files
genes.xlsx <- snakemake@input[["genes"]]
te.xlsx <- snakemake@input[["te"]]

# create output directories
outdir_genes <- snakemake@output[["genes"]]
outdir_te <- snakemake@output[["te"]]

dir.create(outdir_genes, showWarnings = FALSE)
dir.create(outdir_te, showWarnings = FALSE)

# get plotting paramaters
fdr <- -log10(snakemake@params[["fdr"]])
fc <- snakemake@params[["fc"]]

# function to plot volcano for each sheet (i.e. DESeq2 contrast) in given xlsx
volcano <- function(xlsx, outdir){

  for (i in excel_sheets(xlsx)){
    # load data for contrast
    df <- read_excel(xlsx, sheet = i)
  
    # get contrast name
    contrast <- i
    print(paste0("Generating volcano for ",contrast," in ",outdir,"..."))
    
    # log transform padj
    df$log.padj <- -log10(df$padj)
  
    # remove genes with NA in log.padj
    df <- df[!(is.na(df$log.padj)), ]
    
    # add fill colour based on log2FC and pvalue
    df <- df %>%
      mutate(colour = case_when(
        log2FoldChange > fc & log.padj > fdr ~ "red",
        log2FoldChange < -fc & log.padj > fdr ~ "navy",
        log2FoldChange < fc & log.padj < fdr ~ "grey40",
        log2FoldChange > -fc & log.padj < fdr ~ "grey40",
      )) 
    
    # select top 5 down and up regulated for labels
    df.up <- df %>%
      filter(log2FoldChange > fc) %>%
      filter(log.padj > fdr)
    if (nrow(df.up) > 5){
      df.up <- df.up %>% slice(1:5)
    }
      
    df.down <- df %>%
      filter(log2FoldChange < -fc) %>%
      filter(log.padj > fdr)
      if (nrow(df.down) > 5){
        df.down <- df.down %>% slice(1:5)
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
      geom_vline(xintercept = c(-fc,fc),
                 linetype = "dashed", 
                 color = "red", 
                 linewidth = 0.5) +
      geom_hline(yintercept = fdr,
                 linetype = "dashed", 
                 color = "red", 
                 linewidth = 0.5) +
      ggtitle(contrast) 
    
    if (str_detect(xlsx,"deseq2_te")){
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
      
      # save plot to file
      ggsave(file.path(outdir,paste0(contrast,".pdf")), p)
  }
}

# run plotting function for genes and te
volcano(genes.xlsx, outdir_genes)
volcano(te.xlsx, outdir_te)


sink(log, type = "output")
sink(log, type = "message")


