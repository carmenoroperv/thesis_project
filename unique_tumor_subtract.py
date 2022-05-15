import os
import subprocess
from subprocess import Popen
import sys

others_list = snakemake.params.others_list
pt = snakemake.params.pt
output_f = snakemake.params.output_f

input1 = "data/" + str(pt) + "/" + str(output_f) + "/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers"
input2 = "data/" + str(others_list[0][0]) + "/" + str(others_list[0][1]) + "/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers"
output1 = str(others_list[0][2])

print(input1)
print(input2)
print(output1)
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print(os.system(f"kmc_tools simple {input1} -ci5 -cx1000000000 {input2} -ci0 -cx1000000000 kmers_subtract {output1} -ci0 -cx1000000000 -cs1000000000"))
sys.stdout.flush()
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print("First substraction done!")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")

for i in range(0, len(others_list)-1):
    print("Next substraction:")
    print(f"other patient: {others_list[i+1][0]}")
    print(f"output file: {others_list[i+1][2]}")
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    sys.stdout.flush()
    
    input1 = str(others_list[i][2])
    input2 = "data/" + str(others_list[i+1][0]) + "/" + str(others_list[i+1][1]) + "/all_germline_filtered_bams_tumor_ci5_cs1e9/unique_tumor_kmers"
    output1 = str(others_list[i+1][2])
    
    print(os.system(f"kmc_tools simple {input1} -ci5 -cx1000000000 {input2} -ci0 -cx1000000000 kmers_subtract {output1} -ci0 -cx1000000000 -cs1000000000"))
    sys.stdout.flush()
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
   
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")

