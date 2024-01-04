import pandas as pd
import gzip
import shutil
import os
import zipfile
import sys


sys.stderr = open(snakemake.log[0], "w")

# Renaming the qiime2 artifacts for viewing to zip files so the information can be accessed

# Iterating over the input files and renaming them to zip
i = 0
zipped_files = []
while i < len(snakemake.input):
    file = snakemake.input[i]
    name = os.path.splitext(file)[0]
    shutil.copy(file, name + ".zip")
    zipped_files.append(name + ".zip")
    i = i + 1
z = 0
# opening the zip files and saving them to the output directory
while z < len(zipped_files):
    with zipfile.ZipFile(zipped_files[z], "r") as zip_ref:
        name = zipped_files[z].split("/")[-1]
        new_dir = str(snakemake.output) + "/" + name
        if os.path.exists(new_dir):
            shutil.rmtree(new_dir)
        zip_ref.extractall(os.path.splitext(new_dir)[0] + "/")
        os.remove(zipped_files[z])
    z = z + 1
directory = os.listdir(str(snakemake.output))
# Iterating through the directory holding the unzipped files, moving the file content one folder up.
# Folder created while unzipping can be removed, because it gives no further information
j = 0
while j < len(directory):
    subdir = os.listdir(str(snakemake.output) + "/" + directory[j])
    orig_dir = str(snakemake.output) + "/" + directory[j] + "/" + subdir[0]
    new_dir = str(snakemake.output) + "/" + directory[j]
    for f in os.listdir(orig_dir):
        path = orig_dir + "/" + f
        shutil.move(path, new_dir)
    os.rmdir(orig_dir)
    j = j + 1
