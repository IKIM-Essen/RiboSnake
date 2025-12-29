import os
import pandas as pd
from pathlib import Path
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import move, copy2
from datetime import date, datetime
import sys


sys.stderr = open(snakemake.log[0], "w")

# Creating a metadata sample-sheet, that holds additional information concerning the samples.
# Metadata information is needed to be able to create plots and for metadata-specific analysis.

config = snakemake.config

data_root = Path(os.path.abspath(config["data"]))
in_path = Path(os.path.abspath(config["input"]))
metadata = pd.read_csv(str(snakemake.input), header=0, delimiter=",")
metadata.dropna(axis=1,inplace=True)

# 1. collect FASTQ files
incoming = [
    f for f in in_path.glob("*.fastq.gz")
    if f.stat().st_size > 100 and "undetermined" not in f.name.lower()
]

# 2. read date
date = metadata.loc[1, "run_date"]
data_path = data_root / date
data_path.mkdir(parents=True,exist_ok=True)

# 3. copy new files
existing = {f.name for f in data_path.iterdir() if f.is_file()}
to_copy = [f for f in incoming]

for f in to_copy:
    print(f"copying {f} to {data_path}")
    copy2(f, data_path)

# 4. sample names from filenames
sample_list = {f.name.split("_")[0] for f in incoming}

# 5. validate against metadata
metadata_samples = set(metadata["sample_name"])

sample_list = sample_list & metadata_samples

# 6. clean metadata
metadata = metadata.rename(columns={"sample_name": "sample-ID"}).fillna(0)
metadata.to_csv(snakemake.output.metadata, sep="\t", index=False)


# 7. build path assignments
paired = "PairedEnd" in str(config["datatype"])
fastq_map = {s: {"R1": None, "R2": None} for s in sample_list}

for f in to_copy:
    sid = f.name.split("_")[0]
    if sid in fastq_map:
        if "R1" in f.name:
            fastq_map[sid]["R1"] = str(f)
        if "R2" in f.name:
            fastq_map[sid]["R2"] = str(f)

sample_info = metadata.set_index("sample-ID")[["run_date"]].loc[list(sample_list)]

sample_info["path1"] = sample_info.index.map(lambda s: fastq_map[s]["R1"])
if paired:
    sample_info["path2"] = sample_info.index.map(lambda s: fastq_map[s]["R2"])

sample_info.to_csv(snakemake.output.sample_info)