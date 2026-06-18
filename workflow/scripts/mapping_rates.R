# redirect R output to log
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(ggplot2)
library(cowplot)

# get snakemake variables
files <- snakemake@input

# create df for storing mapping rates
df <- as.data.frame(matrix(ncol = 2, nrow = 0))
names(df) <- c("sample", "mapping_rate")

# get mapping rates from STAR log files
counter <- 1

for (x in files) {
  # get sample name
  sample <- system(
    paste0("echo ", x, "| sed 's/Log.final.out//'"),
    intern = TRUE
  )
  sample <- basename(sample)

  # get mapping rate
  rate <- system(
    paste0(
      'grep "Uniquely mapped reads %" ',
      x,
      " | awk '{print $NF}' | sed 's/%//'"
    ),
    intern = TRUE
  )

  # add to df
  df[counter, "sample"] <- sample
  df[counter, "mapping_rate"] <- rate

  counter <- counter + 1
}

# round values to 1 decimal
df$mapping_rate <- as.numeric(df$mapping_rate)
df$mapping_rate <- round(df$mapping_rate, digits = 1)

# create plot
p <- ggplot(df, aes(x = sample, y = mapping_rate)) +
  geom_bar(stat = "identity", fill = "aquamarine4", colour = "black") +
  theme_cowplot(16) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("Mapping rate (%)") +
  xlab(NULL)

# save plot
ggsave(snakemake@output[[1]], p)

# Save mapping rates to CSV
write.csv(df, snakemake@output[["csv"]], row.names = FALSE)

# close redirection of output/messages
sink(log, type = "output")
sink(log, type = "message")
