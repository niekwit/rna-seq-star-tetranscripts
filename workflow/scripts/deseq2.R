# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


library(DESeq2)
library(dplyr)
library(stringr)
library(biomaRt)
library(openxlsx)

# load Snakemake variables
count.files <- snakemake@input[["counts"]]
genome <- snakemake@params[["genome"]]

# Use first count file as template for count matrix
countMatrix <- read.delim(count.files[1]) 
names(countMatrix) <- c("index", sub(".cntTable","",basename(count.files[1])))

# add all count data to countMatrix
for (i in seq(from=2,to=length(count.files))) {
    sample <- sub(".cntTable","",basename(count.files[i]))
    df <- read.delim(count.files[i])
    colnames(df) <- c("index", sample)
    
    countMatrix <- full_join(countMatrix, 
                             df, 
                             by="index")
  }

# remove lines with all 0s
countMatrix <- countMatrix[rowSums(countMatrix[,2:ncol(countMatrix)]) > 0,]

# create named index
rownames(countMatrix) <- countMatrix$index
countMatrix$index <- NULL

#### Load experiment information ####
samples <- read.csv("config/samples.csv", header=TRUE)
genotypes <- unique(samples$genotype)
treatments <- unique(samples$treatment)

if (length(treatments) > 1){
  samples$comb <- paste0(samples$genotype,"_",samples$treatment)
} else {
  samples$comb <- paste0(samples$genotype)
}


# create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = countMatrix,
                              colData = samples,
                              design = ~ comb)

# save DESeqDataSet to file (input for other scripts)
save(dds, file=snakemake@output[["rdata"]])

# load reference samples
references <- unique(samples[samples$reference == "yes" ,]$comb)
if (length(references) == 0){
  stop("ERROR: No reference samples found. Please check your samples.csv file.")
}

# create nested lists to store all pairwise comparisons (top level:references, lower level: samples without reference)
df.list.genes <- vector(mode="list", length=length(references))
for (i in seq_along(references)){
  df.list.genes[[i]] <- vector(mode="list", length=(length(unique(samples$comb)) - 1 ))
}

df.list.te <- vector(mode="list", length=length(references))
for (i in seq_along(references)){
  df.list.te[[i]] <- vector(mode="list", length=(length(unique(samples$comb)) - 1 ))
}

tryCatch({
  # Annotate non-TE Ensembl gene IDs with biomaRt
  genes <- row.names(countMatrix)
  if (grepl("hg",genome, fixed=TRUE)) {
    ensembl <- useEnsembl(biomart="genes", 
                          dataset = "hsapiens_gene_ensembl", 
                          mirror = "www")
    genes <- genes[grepl("ENSG[0-9]{11}+", genes, perl=TRUE)]
    
  } else if (grepl("mm",genome, fixed=TRUE)){
    ensembl <- useEnsembl(biomart="genes", 
                          dataset = "mmusculus_gene_ensembl", 
                          mirror = "www")
    genes <- genes[grepl("ENSMUSG[0-9]{11}+", genes, perl=TRUE)]
  } # add other genomes later
  
  # gene annotation
  gene.info <- getBM(filters = "ensembl_gene_id", 
                     attributes = c("ensembl_gene_id", 
                                    "external_gene_name",
                                    "description",
                                    "gene_biotype", 
                                    "chromosome_name",
                                    "start_position",
                                    "end_position", 
                                    "percentage_gene_gc_content"), 
                     values = genes,
                     mart = ensembl)
}, error = function(error) {
  # biomaRt is very prone to errors, so catch errors here and print a more informative error message
  message("bioMart error. There might be an issue with the servers. Please run Snakemake again later.")
  message("Here's the original error message:")
  message(conditionMessage(error))
  }
)

# performs pair-wise comparisons for each reference sample
for (r in seq(references)){
  # copy dds
  dds_relevel <- dds
  
  # relevel dds to reference sample
  dds_relevel$comb <- relevel(dds$comb, ref = references[r])
  
  # differential expression analysis
  dds_relevel <- DESeq(dds_relevel)
  
  # get all pairwise comparisons
  comparisons <- resultsNames(dds_relevel)
  comparisons <- strsplit(comparisons," ")
  comparisons[1] <- NULL
  
  # create df for each comparison
  for (c in seq(comparisons)){
    
    # get name of comparison
    comparison <- comparisons[[c]]
    comparison <- str_replace(comparison,"comb_","") 
    
    res <- results(dds_relevel, name=comparisons[[c]])
    
    df <- as.data.frame(res) %>%
      mutate(ensembl_gene_id = res@rownames, .before=1)
    
    # get non-TE genes
    if (grepl("hg",genome) ==TRUE) {
      df.genes <- df[grepl("ENSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    } else if (grepl("mm",genome) == TRUE){
      df.genes <- df[grepl("ENSMUSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    } # add other genomes later
    
    # get TE genes
    if (grepl("hg",genome) ==TRUE) {
      df.te <- df[!grepl("ENSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    } else if (grepl("mm",genome) == TRUE){
      df.te <- df[!grepl("ENSMUSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    } # add other genomes later
    
    # add gene annotation to df.genes
    df.genes <- left_join(df.genes,gene.info,by="ensembl_gene_id")
    
    # get normalised read counts for each sample  
    temp <- as.data.frame(counts(dds_relevel, normalized=TRUE))
    temp$ensembl_gene_id <- row.names(temp)
    temp$ensembl_gene_id <- gsub("\\.[0-9]*","",temp$ensembl_gene_id) #tidy up gene IDs
    names(temp)[1:length(dds_relevel@colData@listData$sample)] <- dds_relevel@colData@listData$sample
    
    # add normalised read counts to df.genes
    df.genes <- left_join(df.genes,temp,by="ensembl_gene_id")
    df.genes <- df.genes %>%
      mutate(contrast_name = comparison, .before=1)
    
    # add df.genes to df.list.genes
    df.list.genes[[r]][[c]] <- df.genes
    
    # add normalised read counts to df.te
    df.te <- left_join(df.te,temp,by="ensembl_gene_id")
    df.te <- df.te %>%
      mutate(contrast_name = comparison, .before=1)
    
    # add df.te to df.list.te
    df.list.te[[r]][[c]] <- df.te
    
  }
}

# function to flatten nested lists (https://stackoverflow.com/questions/16300344/how-to-flatten-a-list-of-lists/41882883#41882883)
flattenlist <- function(x){  
  morelists <- sapply(x, function(xprime) class(xprime)[1]=="list")
  out <- c(x[!morelists], unlist(x[morelists], recursive=FALSE))
  if(sum(morelists)){ 
    Recall(out)
  }else{
    return(out)
  }
}

# flatten lists
df.list.genes <- flattenlist(df.list.genes)
df.list.te <- flattenlist(df.list.te)

# get contrast names from each df in list
names.genes <- lapply(df.list.genes, function(x) unique(x$contrast_name))
names(df.list.genes) <- names.genes

names.te <- lapply(df.list.te, function(x) unique(x$contrast_name))
names(df.list.te) <- names.te

# write each df also to separate csv file
save2csv <- function(df.list, type){
  for (i in seq(df.list)){
    write.csv(df.list[[i]], 
              paste0("results/deseq2/", names(df.list)[i], type, ".csv"), 
              row.names = FALSE)
  }
}
save2csv(df.list.genes, "_genes")
save2csv(df.list.te, "_te")

# write output to file
write.xlsx(df.list.genes, 
           snakemake@output[["genes"]],
           colNames = TRUE)

write.xlsx(df.list.te,
           snakemake@output[["te"]],
           colNames = TRUE)


sink(log, type = "output")
sink(log, type = "message")

