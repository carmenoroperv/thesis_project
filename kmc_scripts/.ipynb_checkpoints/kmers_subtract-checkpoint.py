import subprocess


print("Subtracting kmers")

input_tumor_kmers = snakemake.input.tumor_pre
input_tumor_kmers = input_tumor_kmers.split(".")[0]
print(f"Input tumor kmers: {input_tumor_kmers}")

input_germline_kmers = snakemake.input.germline_pre
input_germline_kmers = input_germline_kmers.split(".")[0]
print(f"Input germline kmers: {input_germline_kmers}")

output_file_pre = snakemake.output.unique_tumor_kmers_pre
output_file = output_file_pre.split(".")[0]
print(f"Output files path: {output_file}")

def execute_cmd(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(process.stdout.readline, ""):
        yield stdout_line 
    process.stdout.close()
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

command = f"kmc_tools simple ../{input_tumor_kmers} ../{input_germline_kmers} kmers_subtract ../{output_file}"

print(f"KMC command to execute: {command}")

for path in execute_cmd(command):
    print(path, end="")