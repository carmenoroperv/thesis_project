library(tidyverse)


create_fasta <- function(dataset, pt, folder) {
    file.create(as.character(snakemake@output[["fasta"]]))
    cat(paste0(">", "kmers"), file=as.character(snakemake@output[["fasta"]]), append=TRUE, sep="\n")
    for (row in 1:nrow(dataset)){
        cat(dataset[row, "V10"], file=snakemake@output[["fasta"]], append=TRUE, sep="\n")
    }
    
}

pt = as.character(snakemake@params[["pt"]])
folder = as.character(snakemake@params[["folder"]])
dataset = read.delim(snakemake@input[["tsv_file"]], sep="\t", header=TRUE)
dataset$V10 <- as.character(dataset$V10)

print(pt)
print(head(dataset))
print(dim(dataset))

create_fasta(dataset, pt, folder)

