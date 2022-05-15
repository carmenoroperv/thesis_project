import pandas as pd
import sys

df = pd.read_csv(snakemake.input.kmer_dump, sep='\t', header =  None)

df = df.rename(columns={0: 'kmer', 1: 'count'})
print(df.head())
print(df.shape)

for index, row in df.iterrows():
    header = ">index" + str(index)
    
    if index % 10000 == 0:
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
            
    