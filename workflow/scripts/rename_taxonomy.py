import pandas as pd
import os
import sys


sys.stderr = open(snakemake.log[0], "w")

# Renaming the headers in the taxonomy.tsv file, so that they match with the headers of the table-biom file.
# When the headers are equal, the files can be merged into one file.

taxonomy_file = str(snakemake.input) + "/" + "taxonomy.tsv"
print(os.listdir(str(snakemake.input))[0])
taxonomy_df = pd.read_csv(taxonomy_file, delimiter="\t", header=0)
taxonomy_df.rename(
    columns={"Feature ID": "#OTU ID", "Taxon": "taxonomy", "Consensus": "confidence"},
    inplace=True,
)
taxonomy_df.to_csv(str(snakemake.output), sep="\t", index=False)
