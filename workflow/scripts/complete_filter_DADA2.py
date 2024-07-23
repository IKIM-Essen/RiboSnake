import os
import pandas as pd
import shutil
import zipfile
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import copy, copy2
from datetime import date, datetime
import sys

sys.stderr = open(snakemake.log[0], "w")


dada2 = pd.read_csv(
    str(snakemake.input.dada2) + "/dada2-stats-visual/data/metadata.tsv",
    sep="\t",
    header=0,
    index_col=0,
)
dada2.drop(["#q2:types"], axis=0, inplace=True)

if "percentage of input merged" in dada2.columns:
    dada2.drop(
        [
            "percentage of input passed filter",
            "percentage of input merged",
            "percentage of input non-chimeric",
        ],
        axis=1,
        inplace=True,
    )
else:
    dada2.drop(
        [
            "percentage of input passed filter",
            "percentage of input non-chimeric",
        ],
        axis=1,
        inplace=True,
    )
length = pd.read_csv(
    str(snakemake.input.length)
    + "/table-cluster-lengthfilter/data/sample-frequency-detail.csv",
    sep=",",
    header=0,
    index_col=0,
)
length.rename(columns={"0": "Reads after length filter"}, inplace=True)
before_abundance = pd.read_csv(
    str(snakemake.input.before_abundance) + "/sample-frequency-detail.csv",
    sep=",",
    header=0,
    index_col=0,
)
before_abundance.rename(columns={"0": "Reads before abundance filter"}, inplace=True)
complete = pd.read_csv(
    str(snakemake.input.final) + "/sample-frequency-detail.csv",
    sep=",",
    header=0,
    index_col=0,
)
complete.rename(columns={"0": "Reads after abundance filter"}, inplace=True)

merged_df = pd.concat([dada2, length, before_abundance, complete], axis=1)

merged_df = merged_df.fillna(0)

# Convert all numbers to integers
merged_df = merged_df.apply(pd.to_numeric, errors="ignore", downcast="integer")

# Generate HTML table
html_table = merged_df.to_html(index=True, classes="qiime2-table")

# Add CSS styling
html_content = f"""
<!DOCTYPE html>
<html>
<head>
<style>
.qiime2-table {{
    font-family: Arial, sans-serif;
    border-collapse: collapse;
    width: 100%;
}}
.qiime2-table th, .qiime2-table td {{
    border: 1px solid #dddddd;
    text-align: left;
    padding: 8px;
}}
.qiime2-table th {{
    background-color: #f2f2f2;
}}
</style>
</head>
<body>

{html_table}

</body>
</html>
"""
# Save the HTML content to a file
with open(str(snakemake.output), "w") as f:
    f.write(html_content)
