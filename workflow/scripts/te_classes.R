# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(reshape2)

# Get cut offs
fdr <- as.numeric(snakemake@params["fdr"])
lfc <- as.numeric(snakemake@params["lfc"])

# Load data
files <- snakemake@input[["te_csv"]]

# Empty list to store data
te <- list()

# Count all TE classes for each file for plotting
for (file in files) {
  # read in data
  data <- read.csv(file)
  
  # Extract sample name
  sample <- data %>%
    dplyr::select(contrast_name) %>%
    unique() %>%
    pull()
  
  # Annotate data with TE class
  data <- data %>%
    mutate(class = str_split(ensembl_gene_id, ":", simplify = TRUE)[, 3]) %>%
    mutate(class = case_when(class == "DNA" ~ "Transposon",
                             class == "DNA?" ~ "Transposon?",
                             class == "LTR" ~ "ERV",
                             class == "RC" ~ "Rolling circle",
                             TRUE ~ class)) %>%
    mutate(effect = case_when(log2FoldChange > lfc & padj < fdr  ~ "Upregulated",
                              log2FoldChange < -lfc & padj < fdr ~ "Downregulated")) %>%
    dplyr::filter(effect %in% c("Upregulated", "Downregulated"))
  
  # Count number of genes in each TE class
  df.up <- data %>%
    dplyr::filter(effect == "Upregulated") %>%
    group_by(class) %>%
    summarise(Upregulated = n()) %>%
    mutate(sample = sample)

  df.down <- data %>%
    dplyr::filter(effect == "Downregulated") %>%
    group_by(class) %>%
    summarise(Downregulated = n()) %>%
    mutate(sample = sample)
  
  # merge data and replace NA with zero
  if (nrow(df.up) != 0 & nrow(df.down) != 0 | nrow(df.up) != 0 & nrow(df.down) == 0 ) {
    df <- left_join(df.up, df.down, by = join_by(class, sample))
  } else if (nrow(df.up) == 0 & nrow(df.down) != 0) {
    df <- left_join(df.down, df.up, by = join_by(class, sample))
  } else if (nrow(df.up) == 0 & nrow(df.down) == 0) {
    next
  }
  
  df <- df %>% 
    replace(is.na(.), 0) %>%
    melt(id.vars = c("class", "sample"),
         variable.name = "effect",
         value.name = "Count")
  
  # Add data to list
  te[[sample]] <- df
}

# Function to plot TE classes in bar graph
te_classes <- function(df) {
  # Get sample name
  sample <- df %>%
    pull(sample) %>%
    unique()
  
  # Plot data
  p <- ggplot(df, aes(x = class, 
                      y = Count, 
                      fill = effect)) +
    geom_bar(stat = "identity", 
             position = "dodge", 
             colour = "black") +
    scale_fill_manual(values = c("Upregulated" = "red", 
                                 "Downregulated" = "navy"),
                      name = NULL) +
    labs(x = NULL, 
         y = "Number of differential elements") +
    theme_cowplot(16) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          legend.justification = "center") + 
    labs(title = sample) +
    scale_y_continuous(expand = c(0, 0))
  
  # Save plot
  ggsave(paste0("results/plots/te_classes/", 
         sample, 
         "_te_classes.pdf"), 
         plot = p)
}

# Plot data
lapply(te, te_classes)

# Close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
