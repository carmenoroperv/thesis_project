import os

others_list = snakemake.params.others_list
pt = snakemake.params.pt
output_f = snakemake.params.output_f

os.system(f"kmc_tools simple data/{pt}/{output_f}/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers -ci5 -cx1000000000 data/{others_list[0][0]}/{others_list[0][1]}/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers -ci0 -cx1000000000 kmers_subtract {others_list[0][2]} -ci0 -cx1000000000 -cs1000000000") 

print("First substraction done!")
for i in range(0, len(others_list)-1):
    print("Next substraction:")
    print(f"other patient: {others_list[i+1][0]}")
    print(f"output file: {others_list[i+1][2]}")
    os.system(f"kmc_tools simple {others_list[i][2]} -ci5 -cx1000000000 data/{others_list[i+1][0]}/{others_list[i+1][1]}/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers_substraction/unique_tumor_kmers -ci0 -cx1000000000 kmers_subtract {others_list[i+1][2]} -ci0 -cx1000000000 -cs1000000000") 
