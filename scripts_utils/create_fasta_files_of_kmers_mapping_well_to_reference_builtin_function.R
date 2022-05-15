library(tidyverse)
#install.packages("seqinr", repos = "http://cran.us.r-project.org/src/contrib/seqinr_4.2-8.tar.gz")
library(seqinr)

pt = as.character(snakemake@params[["pt"]])
folder = as.character(snakemake@params[["folder"]])
dataset = read.delim(snakemake@input[["tsv_file"]], sep="\t", header=TRUE)
dataset$V10 <- as.character(dataset$V10)
seq <- dataset$V10
names <- dataset$V1

print(pt)
print(head(dataset))
print(dim(dataset))

write.fasta(seq, names, snakemake@output[["fasta"]], open = "w", nbchar = 50, as.string = TRUE)

