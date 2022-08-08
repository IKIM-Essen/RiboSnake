import pandas as pd
import sys


sys.stderr = open(snakemake.log[0], "w")

fasta = open(str(snakemake.input), "r")
fasta_upper_file = open(str(snakemake.output), "w")
print("read!")
for x in fasta.read():
    fasta_upper = x.upper()
    fasta_upper_file.write(fasta_upper)
print("ready")
fasta.close()
fasta_upper_file.close()
