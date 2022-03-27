import subprocess


print("Subtracting kmers")

input_kmers = snakemake.input.kmers_pre
input_kmers = input_kmers.split(".")[0]
print(f"Input kmers: {input_kmers}")

output_file = snakemake.output.plotdata
print(f"Output file: {output_file}")

def execute_cmd(cmd):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(process.stdout.readline, ""):
        yield stdout_line 
    process.stdout.close()
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

command = f"kmc_tools transform ../{input_kmers} histogram ../{output_file}"

print(f"KMC command to execute: {command}")

for path in execute_cmd(command):
    print(path, end="")
