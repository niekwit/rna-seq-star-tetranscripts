# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


library(DESeq2)
library(dplyr)
library(stringr)
library(biomaRt)
library(openxlsx)

# create output dir
dir.create("deseq2", showWarnings=FALSE)

#### load snakemake variables####
count.files <- snakemake@input[["counts"]]
gtf <- snakemake@input[["gtf"]]
all.genes.txt <- snakemake@input[["genes"]]
genome <- snakemake@params[["genome"]]

#### Load gene information####

genes <- scan(all.genes.txt, what="character")


#### Create count matrix for DESeq2 from count files ####
countMatrix <- data.frame(matrix(ncol=1, nrow=length(genes)))
names(countMatrix) <- c("index")
countMatrix$index <- genes
rownames(countMatrix) <- countMatrix$index

# add count data to countMatrix
for (i in count.files) {
    sample <- basename(i)
    sample <- sub(".cntTable","",sample)
    df <- read.delim(i)
    colnames(df) <- c("index", sample)
    
    #get TE gene name
    df$index <- str_split_fixed(df$index,":", n=2)[,1]
    
    countMatrix <- left_join(countMatrix, df, by="index")
  
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

# create nested lists to store all pairwise comparisons (top level:references, lower level: samples without reference)
df.list.genes <- vector(mode="list", length=length(references))
for (i in seq_along(references)){
  df.list.genes[[i]] <- vector(mode="list", length=(length(unique(samples$comb)) - 1 ))
}

df.list.te <- vector(mode="list", length=length(references))
for (i in seq_along(references)){
  df.list.te[[i]] <- vector(mode="list", length=(length(unique(samples$comb)) - 1 ))
}

# load data for gene annotation
mart <- useMart("ensembl")

if (grepl("hg",genome, fixed=TRUE)) {
  mart <- useDataset("hsapiens_gene_ensembl", mart = mart)
} else if (grepl("mm",genome, fixed=TRUE)){
  mart <- useDataset("mmusculus_gene_ensembl", mart = mart)
} # add other genomes later

for (r in seq(references)){
  
  # copy dds
  dds_relevel <- dds
  
  # for each reference sample, set it as reference in dds
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
    df.genes <- df[grepl("ENSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    
    # get TE genes
    df.te <- df[!grepl("ENSG[0-9]{11}+",df$ensembl_gene_id, perl=TRUE),] 
    
    # get gene annotation
    gene.info <- getBM(filters = "ensembl_gene_id", 
                       attributes = c("ensembl_gene_id", 
                                      "external_gene_name",
                                      "description",
                                      "gene_biotype", 
                                      "chromosome_name",
                                      "start_position",
                                      "end_position", 
                                      "percentage_gene_gc_content"), 
                       values = df.genes$ensembl_gene_id, 
                       mart = mart)
    
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

# write output to file
write.xlsx(df.list.genes, 
           snakemake@output[["genes"]],
           colNames = TRUE)

write.xlsx(df.list.te,
           snakemake@output[["te"]],
           colNames = TRUE)


sink(log, type = "output")
sink(log, type = "message")



