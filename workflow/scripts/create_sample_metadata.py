import os
import pandas as pd
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import move, copy2
from datetime import date, datetime
# Creating a metadata sample-sheet, that holds additional information concerning the samples.
# Metadata information is needed to be able to create plots and for metadata-specific analysis.

DATA_PATH = snakemake.input.data
IN_PATH = snakemake.input.incoming
today = date.today().strftime("%Y-%m-%d")
incoming_files = []
# Running through all files in the IN_PATH and checking them
# Adding the files to a list, when they meet the requirements
for f in listdir(IN_PATH):
    if (
        path.isfile(path.join(IN_PATH, f))
        and "ndetermined" not in f
        and ".fastq.gz" in f
        and os.stat(IN_PATH + f).st_size > 100
    ):
        incoming_files.append(f)
    else:
        print(f, "not used")
# Adding a direcory for the specific date as subdirectory to the DATA_PATH
DATA_PATH += snakemake.params.date
if not os.path.isdir(DATA_PATH):
    mkdir(DATA_PATH)
# Checking if some files are already present in the DATA_PATH, those are not copied
data_files = [f for f in listdir(DATA_PATH) if path.isfile(path.join(DATA_PATH, f))]
files_not_to_copy = [f for f in data_files if f in incoming_files]
files_to_copy = [f for f in incoming_files if f not in data_files]
for file in files_to_copy:
    if file.endswith(".fastq.gz") and not "ndetermined" in file:
        move(IN_PATH + file, DATA_PATH)

# Reading the sample names and the metadata from the file-names and the metadata.csv file, that needs to be provided
files = os.listdir(IN_PATH)
print(files)
sample_list = []
for name in files:
    sample = name.split("_")[0]
    sample_list.append(sample)
sample_list = list(set(sample_list))
metadata = pd.read_csv(snakemake.input.metadata, header = 0, delimiter = ",")
# Samples, that are not mentioned in the metadata.csv are excluded from the sample-metadata-sheet
for name in sample_list:
    if name not in metadata["sample_name"]:
        print("No metadata was found for " + name)
        sample_list.remove(name)
        # Also remove the fastq file from folder?
# Replacing empty metadata-columns with "NaN" and after that removing them from the file
# Only columns holding information should be added to the sample-metadata-sheet
nan_value = float("NaN")
metadata.replace("", nan_value, inplace=True)
metadata.dropna(how='any', axis=1, inplace=True)
metadata.rename(columns = {"sample_name":"sample-ID"}, inplace = True)
metadata.to_csv(snakemake.output, sep = '\t', index=False)

