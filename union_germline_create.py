import os
import sys

samples = snakemake.params.samples
folders = snakemake.params.folders

print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print(f"Combining germlines: FIRST COMBINATION")
print(f"patient1: {samples[0]}")
print(f"folder1: {folders[0]}")
print(f"patient2: {samples[1]}")
print(f"folder2: {folders[1]}")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
sys.stdout.flush()

input1 = "data/" + str(samples[0]) + "/" + str(folders[0]) + "/cs1e9/germline_kmers_cs1e9"
input2 = "data/" + str(samples[1]) + "/" + str(folders[1]) + "/cs1e9/germline_kmers_cs1e9"
output1 = "data/union_germline_combining/germline_kmers_cs1e9_" + str(samples[0]) + "__" + str(samples[1])
print(f"Input1: {input1}")
print(f"Input2: {input2}")
print(f"Output1: {output1}")
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
print(os.system(f"kmc_tools simple {input1} -ci0 -cx1000000000 {input2} -ci0 -cx1000000000 union {output1} -ci0 -cx1000000000 -cs1000000000"))
sys.stdout.flush()
print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")


samples_sub = samples[2:]
folders_sub = folders[2:]

for i, (pt, folder) in enumerate(zip(samples_sub, folders_sub)):
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    print(f"Combining germlines: {i}")
    print(f"patient: {pt}")
    print(f"folder: {folder}")
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
    sys.stdout.flush()
    
    if i == 0: 
        input1 = "data/union_germline_combining/germline_kmers_cs1e9_" + str(samples[0]) + "__" + str(samples[1])
    else: 
        input1 = output1
    input2 = "data/" + str(pt) + "/" + str(folder) + "/cs1e9/germline_kmers_cs1e9"
    
    if i == len(samples_sub)-1:
        output1 = "data/germline_union_NEW"
    else:
        output1 = input1 + "__" + str(pt)
    
    print(f"Input1: {input1}")
    print(f"Input2: {input2}")
    print(f"Output1: {output1}")
    
    print(os.system(f"kmc_tools simple {input1} -ci0 -cx1000000000 {input2} -ci0 -cx1000000000 union {output1} -ci0 -cx1000000000 -cs1000000000"))
    sys.stdout.flush()
   
    print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")

