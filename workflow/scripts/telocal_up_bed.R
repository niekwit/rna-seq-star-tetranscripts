# Redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)

fdr <- snakemake@params[["fdr"]]
lfc <- snakemake@params[["fc"]]

de <- read.csv(snakemake@input[["csv"]])

# TElocal locations file: TE name, then chromosome:start-stop:strand (1-based, inclusive)
loc <- read.delim(snakemake@input[["locations"]], check.names = FALSE)
names(loc)[1:2] <- c("TE", "location")
loc <- loc %>%
  extract(
    location,
    into = c("chrom", "start", "end", "strand"),
    regex = "^(.+):([0-9]+)-([0-9]+):([+-])$",
    convert = TRUE
  )

up <- de %>%
  filter(!is.na(padj), log2FoldChange > lfc, padj < fdr)

# TElocal row names are <locus>:<gene>:<family>:<class>; the locations file uses <locus>
up <- up %>% mutate(te = sub(":.*$", "", ensembl_gene_id))
unmatched <- sum(!(up$te %in% loc$TE))
print(paste(unmatched, "upregulated TEs have no location and are skipped"))

bed <- up %>%
  inner_join(loc, by = c("te" = "TE")) %>%
  mutate(
    start = start - 1,
    score = pmin(1000, round(-log10(padj) * 100))
  ) %>%
  select(chrom, start, end, name = ensembl_gene_id, score, strand) %>%
  arrange(chrom, start)

write.table(
  bed,
  snakemake@output[["bed"]],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

sink(log, type = "output")
sink(log, type = "message")
