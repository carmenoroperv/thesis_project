library(ggplot2)
library(tidyverse)


data=read.table(snakemake@input[["input_filename"]])
data_new1 <- data %>% filter(V1 <= 200)


qplot(V1, V2,data=data_new1, geom="col", xlab="K-mer count", ylab = "Count")
ggsave(snakemake@output[["output_filename"]])