# redirect R output to log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggfortify)
library(RColorBrewer)
library(GenomicFeatures)

# create output dir
dir.create("plots-salmon", showWarnings=FALSE)

# load trx to gene info
gtf <- snakemake@params[["gtf"]]

txdb.filename <- paste0(gsub(".gtf","",gtf, fixed=TRUE),".txdb")

if(file.exists(txdb.filename) == FALSE){
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  txdb <- loadDb(txdb.filename)
} else {txdb <- loadDb(txdb.filename)}

# create df with just trx names
k <- keys(txdb,keytype="TXNAME")
df <- AnnotationDbi::select(txdb,k,"GENEID","TXNAME") %>%
  dplyr::select(TXNAME) %>%
  rename(Name = TXNAME) %>%
  as.data.frame()

# load all Salmon quant.sf file and merge to df with name column
files <- snakemake@input[["salmon"]]
for (i in files){
  
  sample_name <- gsub("/quant.sf","",i, fixed=TRUE) %>% 
    gsub("salmon/","",., fixed=TRUE)
  
  # load salmon quant file, select name and TPM columns, rename TPM column to sample name
  tmp <- read.delim(i, header=TRUE) %>%
    dplyr::select(Name, TPM) %>%
    rename(!!sample_name := TPM) 
  
  # join to df
  df <- left_join(df, tmp, by="Name") 
  
}

# remove rows with only zeros
df <- df %>%
  dplyr::select(-c(Name)) %>%
  filter(rowSums(.) != 0) %>%
  t() %>%
  as.data.frame()

# PCA
pca <- prcomp(df, scale=TRUE)

# load sample info
samples <- read.csv("samples.csv", header=TRUE)

samples.number <- length(samples$sample)
genotype.number <- length(unique(samples$genotype))
condition.number <- length(unique(samples$condition))
replicate.number <- samples.number / genotype.number / condition.number  # assumes same replicate number for all samples


df$Genotype <- samples$genotype
samples$geno_cond <- paste0(samples$genotype, "_", samples$condition)
df$replicate <- as.factor(rep(1:replicate.number, length(unique(samples$geno_cond))))

# check if replicate number is integer
if (replicate.number %% 1 != 0){
  
  print("Replicate number is not an integer. Check samples.csv file.")

  } else {
  
  shape.list <- c(21,22,24,25,23)
  shape.list <- c(15,16,17,18,4,3,8,11,10,13,21,22,23,24,25)
  shape.value.vector <- rep(shape.list[1:replicate.number],genotype.number)
  
  if (condition.number == 1) {
    
    p <- autoplot(pca, 
              data = df,
              size = 8,
              colour = "Genotype",
              shape = "replicate"
              ) +
  theme_set(theme_bw()) +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size=12)
  ) +
  scale_shape_manual(name = "Replicate",
                     values = shape.value.vector) +
  scale_colour_manual(name = "Genotype",
                    values = RColorBrewer::brewer.pal(genotype.number, "Paired")) +
  coord_fixed(ratio = 1)

  
  }
  
}

# save plot to file
ggsave(snakemake@output[[1]], p)


# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
