library(tidyverse)
data_folder <- "all_germline_filtered_bams_tumor_ci5_cs1e9"
paramspace <- read.csv("data/metadata/paramspace_cfDNA_phaseI.csv")
paramspace <- paramspace %>% select(pt_id, cfDNA_folder)
patients <- as.character(unique(paramspace$pt_id))
cfDNA_counts <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(cfDNA_counts) <-c("pt_id", "cfDNA_sample", "count_before_sub", "count_after_sub")




for (i in patients){
    paramspace_patient <- paramspace %>% filter(pt_id == i)
    for (j in paramspace_patient$cfDNA_folder){
        after_sub <- read.csv(paste("data/", i, "/", j, "/all_germline_filtered_bams_tumor_ci5_cs1e9/cfDNA_kmers_substracted_tumor_and_union_germline_histogram_filtered.txt", sep = ""), header=FALSE, sep = "\t")
        sum_after_val <- sum(after_sub[2])
        before_sub <- read.csv(paste("data/", i, "/", j, "/cs1e9/plotdata_cfDNA_cs1e9_filtered.txt", sep = ""), header=FALSE, sep = "\t")
        sum_before_val <- sum(before_sub[2])
        cfDNA_counts[nrow(cfDNA_counts)+1,] = c(i, j, sum_before_val, sum_after_val)
        print("Patient:")
        print(i)
        print("ctDNA sample:")
        print(j)
        flush.console()
    }}



cfDNA_counts = cfDNA_counts %>% mutate(ratio = count_after_sub/count_before_sub)
write.csv(cfDNA_counts, paste("data/Batch_effect_check/", data_folder, "/batch_effect_check.csv", sep = ""), row.names = FALSE)
plot <- ggplot(data = cfDNA_counts)+
    geom_point(aes(x = cfDNA_sample, y=ratio))
ggsave(paste("data/Batch_effect_check/", data_folder, "/batch_effect_check.png", sep = ""))

