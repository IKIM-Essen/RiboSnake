import pandas as pd
import gzip
import shutil
import os
import zipfile
import sys


sys.stderr = open(snakemake.log[0], "w")

# Opening the qiime2 artifacts so the information can be accessed. They already ARE
# zip archives, so they are read in place.
#
# They used to be copied to a sibling "<name>.zip" first, and that copy was read and
# deleted afterwards. The temp name was derived solely from the input path, so two
# rules processing the SAME artifact wrote and removed the identical file underneath
# each other: table-cluster-lengthfilter.qzv is the input of both this rule
# (unzip_frequency_length) and abundance_frequency (relative_abundance.py, same copy
# pattern), and neither depends on the other, so Snakemake runs them in parallel.
# Reading a copy the other job was re-writing raised EOFError -- the bigger the
# artifact, the wider the window (<=40 MB always passed, >=78 MB failed
# reproducibly). Reading the artifact directly removes the shared temp file, a full
# copy of I/O and the doubled disk usage.

# Iterating over the input files and extracting each into the output directory
os.makedirs(str(snakemake.output), exist_ok=True)
for file in snakemake.input:
    stem = os.path.splitext(os.path.basename(file))[0]
    new_dir = str(snakemake.output) + "/" + stem
    if os.path.exists(new_dir):
        shutil.rmtree(new_dir)
    with zipfile.ZipFile(file, "r") as zip_ref:
        zip_ref.extractall(new_dir + "/")

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
