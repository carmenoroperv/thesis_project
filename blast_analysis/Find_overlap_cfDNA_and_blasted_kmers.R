library(tidyverse)

correlations <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(correlations) <-c("pt_id", "cfDNA_sample", "no_overlap")

paramspace <- read.csv("data/metadata/paramspace_cfDNA_phaseI_wo4816.csv")
paramspace <- paramspace %>% select(pt_id, cfDNA_folder)
patients <- as.character(unique(paramspace$pt_id))


for (i in patients){
    print(i)
    flush.console()
    paramspace_cfDNA <- paramspace %>% filter(pt_id == i)   
    for (j in paramspace_cfDNA$cfDNA_folder){
        print(j)
        flush.console()
        unique_cfDNA_table <- read.table(paste("data/", i, "/", j, "/all_germline_filtered_bams_tumor_ci5_cs1e9/cfDNA_kmers_unique_tumor_kmers_intersect.txt", sep = ""), sep="")
        print(head(unique_cfDNA_table))
        colnames(unique_cfDNA_table) <- c("kmer", "count")

        blasted_kmers_table <- read.table(paste("blasting/results/all_germline_filtered_bams_tumor_ci5_cs1e9/correct_kmers/blasted_kmers_correct_", i, ".csv", sep = ""), sep=",", header = TRUE)
        print(head(blasted_kmers_table))
        blasted_kmers_table_k <- blasted_kmers_table %>% select(ssciname, qseq)
        blasted_kmers_table_k <- blasted_kmers_table_k %>% rename(kmer = qseq)
        print(head(blasted_kmers_table_k))
        flush.console()

        intersection <- inner_join(blasted_kmers_table_k, unique_cfDNA_table, by = "kmer")
        print(head(intersection))
        flush.console()
        
        intersection_organisms_order_species_kmers_matches <- intersection %>%
                          group_by(ssciname) %>%
                          summarise(total_matches = n(), n_kmers = n_distinct(kmer)) %>%
                          mutate(Freq = n_kmers/sum(n_kmers)) %>% arrange(desc(n_kmers))
        
        write.csv(intersection_organisms_order_species_kmers_matches, paste("data/", i, "/", j, "/all_germline_filtered_bams_tumor_ci5_cs1e9/cfDNA_kmers_from_other_organisms_NEWY.csv", sep = ""), row.names = FALSE)
        
        intersection_organisms_order_species_kmers_matches$ssciname <- as.character(intersection_organisms_order_species_kmers_matches$ssciname)
        blast_organisms_CRC_pt <- intersection_organisms_order_species_kmers_matches %>% filter(grepl("Bacteroides fragilis|Escherichia coli|Fusobacterium nucleatum|Gemella morbillorum|Parvimonas micra|Peptostreptococcus stomatis|Streptococcus gallolyticus|Bacteroides dorei|Bosea massiliensis|Bacteroides ovatus|Bacteroides vulgatus|Human papillomavirus|Merkel cell polyomavirus|Cytomegalovirus", ssciname))
        write.csv(blast_organisms_CRC_pt, paste("data/", i, "/", j, "/all_germline_filtered_bams_tumor_ci5_cs1e9/cfDNA_kmers_from_crc_organisms_NEWY.csv", sep = ""), row.names = FALSE)
        
        no_overlap <- nrow(intersection)
        correlations[nrow(correlations)+1,] = c(i, j, no_overlap)
    }
}
                         


print(str(correlations))
correlations$no_overlap <- as.numeric(correlations$no_overlap)
write.csv(correlations, snakemake@output[["overlap"]], row.names = FALSE)  