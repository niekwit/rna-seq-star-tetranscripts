# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
#sink(log, type = "output")
#sink(log, type = "message")

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
files <- snakemake@input[["te_csv"]] # Input
output <- snakemake@output[["pdf"]] # Output

# Empty list to store data
te <- list()

# Function to check for emmpty data
no_data <- function(df) {
  if (nrow(df) == 0) {
    print(paste0("No differential TEs found for ", sample, "..."))
    ggsave(output[grepl(sample, output)], plot = ggplot() + theme_void())
  }
}

# Count all TE classes for each file for plotting
for (i in seq_along(files)) {
  # read in data
  print(paste0("Loading data from ", files[[i]],"..."))
  data <- read.csv(files[[i]])
  
  # DESeq2 data might be empty?
  if (no_data(data) == FALSE) {
    next
  }

  # Extract comparison name from file name
  sample <- str_replace(basename(files[[i]]), "_te.csv", "")
  print("Annotating data...")
  
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

  # Check if data has no lines (no differential TEs)
  # If so just output a message and output
  # an empty PDF file (Snakemake expects output)
  if (no_data(data) == FALSE) {
    next
  }

  # Count number of genes in each TE class
  print("Counting number of genes in each TE class...")
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

  # Merge data and replace NA with zero
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
  te[[i]] <- df
}

# Function to plot TE classes in bar graph
te_classes <- function(df) {
  # Get comparison name
  sample <- df %>%
    pull(sample) %>%
    unique()

  # Plot data
  print(paste0("Plotting TE classes for ", sample, "..."))
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
  
  # Get output file that matches sample value from output
  output_file <- output[grepl(sample, output)]

  # Save plot
  print(paste0("Saving plot to ", output_file, "..."))
  ggsave(output_file, plot = p)
}

# Plot data for each comparison
lapply(te, te_classes)

# Close redirection of output/messages
#sink(log, type = "output")
#sink(log, type = "message")
