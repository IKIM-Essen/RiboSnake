import pandas as pd
import sys


sys.stderr = open(snakemake.log[0], "w")

fasta = open(str(snakemake.input), "r")
fasta_upper_file = open(str(snakemake.output), "w")
for x in fasta.read():
    fasta_upper = x.upper()
    fasta_upper_file.write(fasta_upper)
fasta.close()
fasta_upper_file.close()
