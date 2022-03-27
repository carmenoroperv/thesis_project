intersection = read.table(snakemake@input[["intersection"]])
all_cfDNA = read.table(snakemake@input[["all_cfDNA"]])
unique_tumor = read.table(snakemake@input[["unique_tumor"]])
all_tumor = read.table(snakemake@input[["all_tumor"]])

all_tumor = all_tumor[-c(1, 2), ]
print(head(all_tumor))

intersection_sum = sum(intersection$V2*intersection$V1)
all_cfDNA_sum = sum(all_cfDNA$V2*all_cfDNA$V1)
unique_tumor_sum = sum(unique_tumor$V2*unique_tumor$V1)
all_tumor_sum = sum(all_tumor$V2*all_tumor$V1)

print("intersection_sum")
print(intersection_sum)

print("all_cfDNA_sum")
print(all_cfDNA_sum)

print("unique_tumor_sum")
print(unique_tumor_sum)

print("all_tumor_sum")
print(all_tumor_sum)

ratio = (intersection_sum/all_cfDNA_sum)/(unique_tumor_sum/all_tumor_sum)

write.csv(ratio, snakemake@output[["ratio"]], row.names = FALSE)