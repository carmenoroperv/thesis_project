library(tidyverse)






correlations <- data.frame(matrix(ncol = 1, nrow = 0))



    
unique_cfDNA_table <- read.table(snakemake@input[["unique_cfDNA"]], sep="")
blasted_kmers_table <- read.table(snakemake@input[["blasted_kmers"]], sep=",")

blasted_kmers_table_k <- as.data.frame(blasted_kmers_table[-1,])

colnames(unique_cfDNA_table) <- c("kmer", "no")
unique_cfDNA_table <- unique_cfDNA_table %>% select(kmer)

colnames(blasted_kmers_table_k) <- c("kmer")

intersection <- inner_join(blasted_kmers_table_k, unique_cfDNA_table, by = "kmer")

no_overlap <- nrow(intersection)

correlations[nrow(correlations)+1,] = no_overlap
                         


write.csv(correlations, snakemake@output[["no_output"]], row.names = FALSE)    