library(tidyverse)

bag_intersection = read.table(snakemake@input[["bag_intersection"]])
bag_sum = read.table(snakemake@input[["bag_sum"]])

bag_intersection$V2 <- as.numeric(bag_intersection$V2)
bag_sum$V2 <- as.numeric(bag_sum$V2)

bag_intersection_sum = sum(bag_intersection$V2*bag_intersection$V1)
bag_intersection_sum = as.numeric(bag_intersection_sum)

bag_sum_sum = sum(bag_sum$V2*bag_sum$V1)
bag_sum_sum = as.numeric(bag_sum_sum)


print("bag_intersection_sum")
str(bag_intersection_sum)
print(bag_intersection_sum)

print("bag_sum_sum")
str(bag_sum_sum)
print(bag_sum_sum)

ratio = bag_intersection_sum/bag_sum_sum

print("ratio")
print(ratio)

# CI 
bag_intersection_sim = (rpois(10000, bag_intersection_sum))
bag_sum_sim = (rpois(10000, bag_sum_sum))


ratios_sim = bag_intersection_sim/bag_sum_sim
ratios_sim = sort(ratios_sim)

print(length(ratios_sim))
print(head(ratios_sim))
print(tail(ratios_sim))

lowerCI = quantile(ratios_sim,.025)
upperCI = quantile(ratios_sim,.975)

res = tibble(ratio = ratio, lower_CI = lowerCI, upper_CI = upperCI, intersection_sum = bag_intersection_sum, bag_sum_sum = bag_sum_sum)

write.csv(res, snakemake@output[["ratio"]], row.names = FALSE)