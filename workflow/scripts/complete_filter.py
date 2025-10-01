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

config = snakemake.config
datatype = str(config["datatype"])

samples = pd.read_csv(
    str(snakemake.input.samples) + "/paired-seqs/data/per-sample-fastq-counts.tsv",
    sep="\t",
    header=0,
    index_col=0,
)
if datatype == "SampleData[PairedEndSequencesWithQuality]":
    samples.drop(
        ["reverse sequence count"],
        axis=1,
        inplace=True,
    )
samples.rename(columns={"forward sequence count": "Raw reads"}, inplace=True)
trimmed = pd.read_csv(
    str(snakemake.input.trimmed) + "/trimmed-seqs/data/per-sample-fastq-counts.tsv",
    sep="\t",
    header=0,
    index_col=0,
)
if datatype == "SampleData[PairedEndSequencesWithQuality]":
    trimmed.drop(
        ["reverse sequence count"],
        axis=1,
        inplace=True,
    )
trimmed.rename(columns={"forward sequence count": "Reads after trimming"}, inplace=True)
if datatype == "SampleData[PairedEndSequencesWithQuality]":
    joined = pd.read_csv(
        str(snakemake.input.joined) + "/joined-seqs/data/per-sample-fastq-counts.tsv",
        sep="\t",
        header=0,
        index_col=0,
    )
    joined.rename(columns={"forward sequence count": "Reads after joining"}, inplace=True)

first = pd.read_csv(
    str(snakemake.input.first) + "/metadata.tsv", sep="\t", header=0, index_col=0
)
first.drop(["#q2:types"], axis=0, inplace=True)
first.drop(
    [
        "reads-truncated",
        "reads-too-short-after-truncation",
        "reads-exceeding-maximum-ambiguous-bases",
    ],
    axis=1,
    inplace=True,
)
human = pd.read_csv(str(snakemake.input.human), sep=",", header=0, index_col=0)

if "difference" in human.columns:
    human.drop(["difference"], axis=1, inplace=True)

wo_chimera = pd.read_csv(
    str(snakemake.input.wo_chimera)
    + "/table-nonchimeric-wo-borderline/data/sample-frequency-detail.csv",
    sep=",",
    header=0,
    index_col=0,
)
wo_chimera.rename(columns={"0": "Reads after chimera filtering"}, inplace=True)
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

if datatype == "SampleData[PairedEndSequencesWithQuality]":
    merged_df = pd.concat(
        [
            samples,
            trimmed,
            joined,
            first,
            human,
            wo_chimera,
            length,
            before_abundance,
            complete,
        ],
        axis=1,
    )
elif datatype == "SampleData[SequencesWithQuality]":
    merged_df = pd.concat(
        [
            samples,
            trimmed,
            first,
            human,
            wo_chimera,
            length,
            before_abundance,
            complete,
        ],
        axis=1,
    )

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
