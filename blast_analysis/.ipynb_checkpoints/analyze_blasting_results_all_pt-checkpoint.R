library(tidyverse)

pt <- snakemake@params[["patient"]] 
print(paste0("Running blast results analysis for patient: ", pt))
flush.console()


unmapped_kmers_pt <- read.table(snakemake@input[["unmapped_kmers_pt"]], header = FALSE, sep = "")
unmapped_kmers_pt$V1 <- as.character(unmapped_kmers_pt$V1)

indexes = NULL
kmers = NULL
for (row in 1:nrow(unmapped_kmers_pt)){
    val = unmapped_kmers_pt[row, "V1"]
    if (substr(val,1,1) == ">"){
        indexes <- c(indexes, val)
    } else {
        kmers <- c(kmers, val)
    }
}
unmapped_kmers_pt_df <- tibble(index = indexes, kmer = kmers)
n_kmers_unmapped <- length(unique(unmapped_kmers_pt_df$kmer))

print(paste0("There are in total ", n_kmers_unmapped, " k-mers that did not map to the reference"))
flush.console()

blast_pt <- read.table(snakemake@input[["blast_pt"]], header = FALSE, sep = "\t")
colnames(blast_pt) <- c("qseqid", "sseqid", "staxid", "ssciname", "pident", "nident", "length", "mismatch", "gapopen", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")

blast_pt$qseq <- as.character(blast_pt$qseq)
blast_pt$qseqid <- as.character(blast_pt$qseqid)
blast_pt$sseq <- as.character(blast_pt$sseq)
blast_pt$ssciname <- as.character(blast_pt$ssciname)

blast_pt_kmers <- blast_pt %>% mutate(n_char = nchar(qseq))
blast_pt_kmers_correct <- blast_pt_kmers %>% filter(n_char == 50)
write.csv(blast_pt_kmers_correct, snakemake@output[["correct_kmers"]])

n_kmers_blasted <- length(unique(blast_pt_kmers_correct$qseqid))
n_species_matched <- length(unique(blast_pt_kmers_correct$ssciname))

print(paste0("There are in total ", n_kmers_blasted, " k-mers that found at least one match with blast"))
print(paste0("There are in total ", n_species_matched, " species that showed a match with at least one of the k-mers"))
flush.console()
blast_pt_kmers_correct$ssciname <- as.character(blast_pt_kmers_correct$ssciname)
blast_organism_order_pt_species_kmers_matches <- blast_pt_kmers_correct %>%
                          group_by(ssciname) %>%
                          summarise(total_matches = n(), n_kmers = n_distinct(qseq)) %>%
                          mutate(Freq = n_kmers/sum(n_kmers)) %>% arrange(desc(n_kmers))
write.csv(blast_organism_order_pt_species_kmers_matches,  snakemake@output[["all_organisms"]])

print(paste0("The most common organisms found: "))
print(blast_organism_order_pt_species_kmers_matches[1:20,])
flush.console()

blast_organism_order_pt_species_kmers_matches$ssciname <- as.character(blast_organism_order_pt_species_kmers_matches$ssciname)
blast_organisms_CRC_pt <- blast_organism_order_pt_species_kmers_matches %>% filter(grepl("Bacteroides fragilis|Escherichia coli|Fusobacterium nucleatum|Gemella morbillorum|Parvimonas micra|Peptostreptococcus stomatis|Streptococcus gallolyticus|Bacteroides dorei|Bosea massiliensis|Bacteroides ovatus|Bacteroides vulgatus|Human papillomavirus|Merkel cell polyomavirus|Cytomegalovirus", ssciname))

write.csv(blast_organisms_CRC_pt, snakemake@output[["crc_organisms"]])

print(paste0("The most common crc organisms found: "))
print(blast_organisms_CRC_pt[1:20,])
flush.console()