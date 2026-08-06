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

# Reading the sample-table. The qiime2 artifact already IS a zip archive, so it is
# read in place -- it used to be copied to a sibling "<name>.zip" first, which
# collided with unzip_frequency_length reading the same artifact through the same
# temp name (see rename_qzv.py for the full story: parallel rules, EOFError on
# large artifacts).
file = str(snakemake.input)
# Extract the artifact to a folder next to it
dir_name = os.path.dirname(file)
name = os.path.splitext(os.path.basename(file))[0]
new_dir = dir_name + "/" + name
if os.path.isdir(new_dir) and os.path.exists(new_dir):
    shutil.rmtree(new_dir)
with zipfile.ZipFile(file, "r") as zip_ref:
    zip_ref.extractall(new_dir + "/")
directory = new_dir
# Moving the folder inventory one folder up
b = 0
subdir = os.listdir(directory)
while b < len(subdir):
    orig_dir = directory + "/" + subdir[b]
    new_dir = directory
    if os.path.isdir(orig_dir):
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
median_of_sums = frequency["Abundance"].median()

endnumber = median_of_sums * abundance
endnumber = int(endnumber)
if endnumber == 0:
    endnumber += 1
with open(str(snakemake.output.abundance), "w") as f:
    f.write(str(endnumber))
