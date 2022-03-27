import pysam
import sys
import pandas as pd

def find_max_quality(in_file, out_file):

    _BASE_QUAL_CUTOFF = 20

    bam_in = pysam.Samfile(in_file, 'rb')

    reads = 0
    qualities = {}
    for read in bam_in.fetch():
        reads += 1
        if reads % 10000000 == 0:
            print(reads)
            sys.stdout.flush()
        for c in list(read.qual):
            phred_qual = ord(c)-33
            if phred_qual in qualities.keys():
                qualities[phred_qual] += 1
            else: 
                qualities[phred_qual] = 1
    
    
    print(f"Total reads analyzed: {reads}")
    print("Quality results")
    print(qualities)
    
    qualities_list = []
    counts_list = []
    for key, value in qualities.items():
        qualities_list.append(key)
        counts_list.append(value)
    
    res = {'Quality': qualities_list, 'Number_of_bases': counts_list}
    df = pd.DataFrame(res, columns=['Quality','Number_of_bases'])
    df.to_csv(out_file, index=False)
    

            
if __name__ == "__main__":
    in_file = snakemake.input.bam_file
    out_file = snakemake.output.qualities
    find_max_quality(in_file, out_file)
    print("Done!")