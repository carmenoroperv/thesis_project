library(tidyverse)

intersection = read.table(snakemake@input[["intersection"]])
all_cfDNA = read.table(snakemake@input[["all_cfDNA"]])
unique_tumor = read.table(snakemake@input[["unique_tumor"]])
all_tumor = read.table(snakemake@input[["all_tumor"]])

intersection$V2 <- as.numeric(intersection$V2)
all_cfDNA$V2 <- as.numeric(all_cfDNA$V2)
unique_tumor$V2 <- as.numeric(unique_tumor$V2)
all_tumor$V2 <- as.numeric(all_tumor$V2)

intersection_sum = sum(intersection$V2)
intersection_sum = as.numeric(intersection_sum)

all_cfDNA_sum = sum(all_cfDNA$V2)
all_cfDNA_sum = as.numeric(all_cfDNA_sum)

unique_tumor_sum = sum(unique_tumor$V2)
unique_tumor_sum = as.numeric(unique_tumor_sum)

all_tumor_sum = sum(all_tumor$V2)
all_tumor_sum = as.numeric(all_tumor_sum)


print("intersection_sum")
str(intersection_sum)
print(intersection_sum)

print("all_cfDNA_sum")
str(all_cfDNA_sum)
print(all_cfDNA_sum)

print("unique_tumor_sum")
str(unique_tumor_sum)
print(unique_tumor_sum)

print("all_tumor_sum")
str(all_tumor_sum)
print(all_tumor_sum)

ratio = (intersection_sum/all_cfDNA_sum)/(unique_tumor_sum/all_tumor_sum)


intersection_sim = (rpois(10000, intersection_sum))
all_cfDNA_sim = (rpois(10000, all_cfDNA_sum))
unique_tumor_sim = (rpois(10000, unique_tumor_sum))
all_tumor_sim = (rpois(10000, all_tumor_sum))

ratios_sim = (intersection_sim/all_cfDNA_sim)/(unique_tumor_sim/all_tumor_sim)
ratios_sim = sort(ratios_sim)

print(length(ratios_sim))
print(head(ratios_sim))
print(tail(ratios_sim))

lowerCI = quantile(ratios_sim,.025)
upperCI = quantile(ratios_sim,.975)

res = tibble(ratio = ratio, lower_CI = lowerCI, upper_CI = upperCI, intersection_sum = intersection_sum, all_cfDNA_sum = all_cfDNA_sum, unique_tumor_sum = unique_tumor_sum, all_tumor_sum = all_tumor_sum)

write.csv(res, snakemake@output[["ratio"]], row.names = FALSE)