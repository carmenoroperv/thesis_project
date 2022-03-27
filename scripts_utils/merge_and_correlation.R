library(tidyverse)






paramspace <- read.table("data/metadata/paramspace_cfDNA_phaseI.csv", sep=",", header =TRUE)

correlations <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(correlations) <-c("pt_id", "cfDNA_sample", "Correlation")

pt <- snakemake@params[['pt_id']]
 
paramspace_cfDNA <- paramspace %>% filter(pt_id == pt)


    
for (j in paramspace_cfDNA$cfDNA_folder){
        tumor <- read.table(paste("data/", pt, "/", j, "/ci5_cs1e9/intersection_unique_tumor_cfDNAintersection.txt", sep=""))
        colnames(tumor) <- c("kmer", "tumor")
        cfDNA <- read.table(paste("data/", pt, "/", j, "/ci5_cs1e9/cfDNA_kmers_unique_tumor_kmers_intersect_singletons_excluded_tumor_4_cs.txt", sep=""))
        colnames(cfDNA) <- c("kmer", paste("cfDNA", j, sep=""))
        tumor_cfDNA <- full_join(tumor, cfDNA, by = "kmer")
        correlations[nrow(correlations)+1,] = c(pt, j, cor(tumor_cfDNA$tumor, tumor_cfDNA$cfDNA))}
                         


write.csv(correlations, snakemake@output[["correlation_output"]], row.names = FALSE)    