import pysam
import sys
import pandas as pd


def FilterReads(in_file, in_qualities_file, out_file):
    
    qualities = pd.read_csv(in_qualities_file)  
    max_quality = qualities['Quality'].max()

    bam_in = pysam.Samfile(in_file, 'rb')
    bam_out = pysam.Samfile(out_file, 'wb', template=bam_in)

    reads = 0
    passing_reads = 0
    filtered_out = 0
    for read in bam_in.fetch():
        reads += 1
        if reads % 10000000 == 0:
            print(reads)
            sys.stdout.flush()
        
        read_new = ""
        passing_read = True
        for base, qual in zip(read.seq, read.qual):
            if ord(qual)-33 == max_quality:
                read_new = read_new + base
            else: 
                read_new = read_new + "N"
                passing_read = False
    
        read.seq = read_new
        if passing_read == True: 
            passing_reads += 1
        else:
            filtered_out += 1
        
        bam_out.write(read)
        
    
    print(f"Total reads analyzed: {reads}")
    print(f"Reads passing the filtering:  {passing_reads}")
    print(f"Reads where at least one base set to N:  {filtered_out}")
            
            
if __name__ == "__main__":
    in_file = snakemake.input.bam_file
    in_qualities_file = snakemake.input.qualities_file
    out_file = snakemake.output.filtered_bam_file
    FilterReads(in_file, in_qualities_file, out_file)
    print("Done!")