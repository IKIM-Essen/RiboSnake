import pandas as pd
import gzip
import shutil
import os
import zipfile
import sys


sys.stderr = open(snakemake.log[0], "w")

# Extracting the number of total features over all samples from the sample-table.
# Multiplying the number with the relative abundance filtering value to create a threshold for the
# qiime filtering.

# Reading the sample-table, creating a zip file
file = str(snakemake.input)
name = os.path.splitext(file)[0]
shutil.copy(file, name + ".zip")
filename = name + ".zip"
# Extract zip files to folder
with zipfile.ZipFile(filename, "r") as zip_ref:
    name = filename.split("/")[-1]
    dir_name = os.path.dirname(str(snakemake.input))
    new_dir = dir_name + "/" + name
    if os.path.isdir(new_dir) and os.path.exists(new_dir):
        shutil.rmtree(new_dir)
    zip_ref.extractall(os.path.splitext(new_dir)[0] + "/")
name = name.split(".")[0]
directory = os.path.dirname(str(snakemake.input)) + "/" + name
# Moving the folder inventory one folder up
b = 0
subdir = os.listdir(directory)
while b < len(subdir):
    orig_dir = directory + "/" + subdir[b]
    new_dir = directory
    for f in os.listdir(orig_dir):
        path = orig_dir + "/" + f
        shutil.move(path, new_dir)
    b = b + 1
# Read the specific csv holding the information, creating a dataframe, adding up all feature frequencies
datadir = str(snakemake.output.feature_table) + "/"
csv = datadir + "sample-frequency-detail.csv"
frequency = pd.read_csv(csv, header=None, delimiter=",")
frequency.columns = ["Sample", "Abundance"]
# number = frequency["Abundance"].sum()
# Creating the abundance threshold and storing it in an output file
abundance = float(str(snakemake.params))

column_sums = frequency.sum()

# Then, calculate the median value of the column sums
median_of_sums = column_sums.median()

endnumber = median_of_sums * abundance
endnumber = int(endnumber)
if endnumber == 0:
    endnumber += 1
with open(str(snakemake.output.abundance), "w") as f:
    f.write(str(endnumber))
