# Redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(DESeq2)

# Load Snakemake variables
count.files <- snakemake@input[["counts"]]
genome <- snakemake@params[["genome"]]

# Use first count file as template for count matrix
countMatrix <- read.delim(count.files[1]) 
names(countMatrix) <- c("index", sub(".cntTable", "", basename(count.files[1])))

# Add all count data to countMatrix
for (i in seq(from = 2, to = length(count.files))) {
    sample <- sub(".cntTable", "", basename(count.files[i]))
    df <- read.delim(count.files[i])
    colnames(df) <- c("index", sample)
    
    countMatrix <- full_join(countMatrix,
                             df,
                             by = "index")
  }

# Remove lines with all 0s
countMatrix <- countMatrix[rowSums(countMatrix[,2:ncol(countMatrix)]) > 0,]

# Create named index
rownames(countMatrix) <- countMatrix$index
countMatrix$index <- NULL

#### Load experiment information ####
samples <- read.csv("config/samples.csv", header = TRUE)
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)

if (length(treatments) > 1) {
  samples$comb <- paste0(samples$genotype, "_", samples$treatment)
} else {
  samples$comb <- samples$genotype
}

# Check if batch column exists
if ("batch" %in% colnames(samples)) {
  batches <- unique(samples$batch)
  if (length(batches) == 1) {
    batches <- 1
  } else {
    samples$batch <- as.factor(samples$batch)
  }
} else {
  batches <- 1
}

# Check if any samples have been omitted from samples.csv and remove from count.files
# This allows for a re-run of DESeq2 without having to re-run the entire pipeline
all_samples <- samples$sample

if (length((all_samples)) != length(count.files)){
  omitted <- count.files[!(basename(gsub(pattern = "\\.cntTable", "", count.files)) %in% all_samples)]
  print("WARNING: Some samples have been omitted from samples.csv. Resuming without these samples:")
  for (i in seq(omitted)){
    print(basename(gsub(pattern = "\\.cntTable", "", omitted[i])))
   }
  # Remove omitted samples from countMatrix
  countMatrix <- countMatrix[, colnames(countMatrix) %in% all_samples]
}

# Create DESeq2 object
if (length(batches) == 1){
  print("Not including batch factor in DESeq2 design...")
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = samples,
                              design = ~comb)
} else {
  print("Including batch factor in DESeq2 design...")
  dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = samples,
                              design = ~batch + comb)
}

# Save DESeqDataSet to file (input for other scripts)
save(dds, file = snakemake@output[["rdata"]])

# Load reference samples
references <- unique(samples[samples$reference == "yes" , ]$comb)
if (length(references) == 0) {
  stop("ERROR: No reference samples found. Please check your samples.csv file.")
}

# Create nested lists to store all pairwise comparisons (top level:references, lower level: samples without reference)
df.list.genes <- vector(mode="list", length = length(references))
for (i in seq_along(references)) {
  df.list.genes[[i]] <- vector(mode = "list",
                               length = (length(unique(samples$comb)) - 1 ))
}

df.list.te <- vector(mode = "list", length = length(references))
for (i in seq_along(references)) {
  df.list.te[[i]] <- vector(mode = "list",
                            length = (length(unique(samples$comb)) - 1 ))
}

# Get gene IDs
genes <- row.names(countMatrix)
if (grepl("hg", genome, fixed = TRUE)) {
  genes <- genes[grepl("ENSG[0-9]{11}+", genes, perl = TRUE)]
} else if (grepl("mm", genome, fixed=TRUE)) {
  genes <- genes[grepl("ENSMUSG[0-9]{11}+", genes, perl = TRUE)]
} else if (genome == "test") {
  genes <- genes[grepl("ENSG[0-9]{11}+", genes, perl = TRUE)]
}

# Get gene annotation
load(snakemake@input[["edb"]])

# Remove duplicate lines
gene.info <- gene.info[!duplicated(gene.info$ensembl_gene_id), ]

# Performs pair-wise comparisons for each reference sample
for (r in seq(references)) {
  # Copy dds
  dds_relevel <- dds

  # Relevel dds to reference sample
  print(paste("Releveling to reference sample: ", references[r], "..."))
  dds_relevel$comb <- relevel(dds$comb, ref = references[r])

  # Differential expression analysis
  dds_relevel <- DESeq(dds_relevel)

  # Get all pairwise comparisons
  comparisons <- resultsNames(dds_relevel)
  comparisons <- strsplit(comparisons, " ")
  comparisons[1] <- NULL

  # Create df for each comparison
  for (c in seq(comparisons)) {
    # Get name of comparison
    comparison <- comparisons[[c]]
    comparison <- str_replace(comparison, "comb_", "") 

    res <- results(dds_relevel, name = comparisons[[c]])

    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = res@rownames, .before = 1)

    # Get non-TE genes
    if (grepl("hg",genome) ==TRUE) {
      df.genes <- df[grepl("ENSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),] 
    } else if (grepl("mm",genome) == TRUE){
      df.genes <- df[grepl("ENSMUSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),] 
    } else if (genome == "test") {
      df.genes <- df[grepl("ENSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),]
    }

    # Get TE genes
    if (grepl("hg",genome) == TRUE) {
      df.te <- df[!grepl("ENSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),] 
    } else if (grepl("mm",genome) == TRUE){
      df.te <- df[!grepl("ENSMUSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),] 
    } else if (genome == "test") {
      df.te <- df[!grepl("ENSG[0-9]{11}+", df$ensembl_gene_id, perl = TRUE),] 
    }

    # Add gene annotation to df.genes
    df.genes <- left_join(df.genes,gene.info, by = "ensembl_gene_id")

    # Get normalised read counts for each sample  
    temp <- as.data.frame(counts(dds_relevel, normalized = TRUE))
    temp$ensembl_gene_id <- row.names(temp)
    temp$ensembl_gene_id <- gsub("\\.[0-9]*", "", temp$ensembl_gene_id) #tidy up gene IDs
    names(temp)[1:length(dds_relevel@colData@listData$sample)] <- dds_relevel@colData@listData$sample

    # Add normalised read counts to df.genes
    df.genes <- left_join(df.genes,temp, by = "ensembl_gene_id")
    df.genes <- df.genes %>%
      mutate(contrast_name = comparison, .before = 1)

    # Add df.genes to df.list.genes
    df.list.genes[[r]][[c]] <- df.genes

    # Add normalised read counts to df.te
    df.te <- left_join(df.te,temp, by = "ensembl_gene_id")
    df.te <- df.te %>%
      mutate(contrast_name = comparison, .before = 1)

    # Add df.te to df.list.te
    df.list.te[[r]][[c]] <- df.te
  }
}

# Function to flatten nested lists (https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists/41882883#41882883)
flattenlist <- function(x) {  
  morelists <- sapply(x, function(xprime) class(xprime)[1] == "list")
  out <- c(x[!morelists], unlist(x[morelists], recursive = FALSE))
  if(sum(morelists)) { 
    Recall(out)
  } else {
    return(out)
  }
}

# Flatten lists
df.list.genes <- flattenlist(df.list.genes)
df.list.te <- flattenlist(df.list.te)

# Get contrast names from each df in list
names.genes <- lapply(df.list.genes, function(x) unique(x$contrast_name))
names(df.list.genes) <- names.genes

names.te <- lapply(df.list.te, function(x) unique(x$contrast_name))
names(df.list.te) <- names.te

# Write each df to separate csv file
save2csv <- function(df.list, type){
  for (i in seq(df.list)) {
    # Check if df is empty
    stopifnot(nrow(df.list[[i]]) > 0)
    # Write to file
    write.csv(df.list[[i]], 
              paste0("results/deseq2/", names(df.list)[i], type, ".csv"), 
              row.names = FALSE)
  }
}
save2csv(df.list.genes, "_genes")
save2csv(df.list.te, "_te")

# Close log
sink(log, type = "output")
sink(log, type = "message")
