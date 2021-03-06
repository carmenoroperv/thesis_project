library(tidyverse)

intersection = read.table(snakemake@input[["intersection"]])
unique_tumor = read.table(snakemake@input[["unique_tumor"]])

intersection$V2 <- as.numeric(intersection$V2)
unique_tumor$V2 <- as.numeric(unique_tumor$V2)

intersection_sum = sum(intersection$V2*intersection$V1)
unique_tumor_sum = sum(unique_tumor$V2*unique_tumor$V1)

print("Intersection sum")
print(intersection_sum)

print("unique_tumor_sum")
print(unique_tumor_sum)

ratio = (intersection_sum/unique_tumor_sum)

print("ratio")
print(ratio)

# CI for small ratio
print(binom_test <- binom.test(intersection_sum, unique_tumor_sum))
print(binom_test$conf.int)
 
print("lowerCI")
print(binom_test$conf.int[1])
print("upperCI")
print(binom_test$conf.int[2])
#lower_CI = NA, upper_CI = NA
res = tibble(ratio = ratio, lower_CI = binom_test$conf.int[1], upper_CI = binom_test$conf.int[2], intersection_sum_small = intersection_sum, unique_tumor_sum_small = unique_tumor_sum)

write.csv(res, snakemake@output[["ratio"]], row.names = FALSE)