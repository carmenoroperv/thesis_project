import subprocess
import os
import shlex

print(os.listdir("./"))

print("Counting kmers")

input_bam_file = snakemake.input.bam
print(f"Input file path: {input_bam_file}")

input_folder = '/'.join(input_bam_file.split("/")[0:3])
print(input_folder)
print(os.listdir(str(input_folder)))

output_file_pre = snakemake.output.kmers_pre
output_file = output_file_pre.split(".")[0]
print(f"Output files path: {output_file}")

if not os.path.exists(snakemake.params.tmpdir):
    os.mkdir(snakemake.params.tmpdir)
print("TMP directory done")
    



command = f'kmc -k50 -m{snakemake.resources.mem_mb} -sm -t{snakemake.threads} -ci0 -fbam {input_bam_file} {output_file} {snakemake.params.tmpdir}'

comm = shlex.split(command)
print(comm)

print(f"KMC command to execute: {command}")

process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell = True, universal_newlines=True)
while True:
    output = process.stdout.readline()
    print(output.strip())
    # Do something else
    return_code = process.poll()
    if return_code is not None:
        print('RETURN CODE', return_code)
        # Process has finished, read rest of the output 
        for output in process.stdout.readlines():
            print(output.strip())
        break


print("Done")

# in the command end:     # 2> {log}


#https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html