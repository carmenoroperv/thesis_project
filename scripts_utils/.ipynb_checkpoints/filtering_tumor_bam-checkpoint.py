import pysam
import sys


def FilterReads(in_file, out_file):

    def read_ok(read):
        """
        read_ok - reject reads with a low quality (<5) base call
        read - a PySam AlignedRead object
        returns: True if the read is ok
        """
        return all([ord(c)-33 >=  _BASE_QUAL_CUTOFF for c in list(read.qual)])


    _BASE_QUAL_CUTOFF = 20

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
        if read_ok(read):
            passing_reads += 1
            bam_out.write(read)
        else:
            filtered_out += 1
    
    print(f"Total reads analyzed: {reads}")
    print(f"Reads passing the filtering:  {passing_reads}")
    print(f"Reads filtered out:  {filtered_out}")
            
            
if __name__ == "__main__":
    in_file = snakemake.input.bam_file
    out_file = snakemake.output.filtered_bam_file
    FilterReads(in_file, out_file)
    print("Done!")