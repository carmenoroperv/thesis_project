import pandas as pd
import sys


pipeline = snakemake.params.pipeline #"all_germline_ci5_cs1e9"
pt = snakemake.params.patient #"C04689"
folder = snakemake.params.folder #"C299A04689D_cfdna_N295-99"
file_name = snakemake.params.filename "cfDNA_kmers_unique_all_germline_tumor_kmers_intersect_singletons_excluded_tumor_4_cs"

#df = pd.read_csv("../data/C04689/C299A04689D_cfdna_N295-99/all_germline_ci5_cs1e9/cfDNA_kmers_unique_all_germline_tumor_kmers_intersect_singletons_excluded_tumor_4_cs.txt", sep='\t', header =  None)
df = pd.read_csv(snakemake.input.kmer_set, sep='\t', header =  None)

df = df.rename(columns={0: 'kmer', 1: 'count'})
print(df.head())

for index, row in df.iterrows():
    header = ">pipeline " +  pipeline + " | patient " + pt + " | folder " + folder +  " | file " + file_name + " | kmer_index_in_file " + str(index) + " | count_in_file " + str(row["count"]) + " |  "
    
    if index % 100000 == 0:
        print(index)
        sys.stdout.flush()
    
    if index == 0:
        with open(snakemake.output.fasta, "w") as file_object:
            file_object.write(header)
            file_object.write("\n")
            file_object.write(row["kmer"])
            file_object.write("\n")
    else: 
        with open(snakemake.output.fasta, "a") as file_object:
            file_object.write(header)
            file_object.write("\n")
            file_object.write(row["kmer"])
            file_object.write("\n")
            
    