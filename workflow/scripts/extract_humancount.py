import os
import shutil
import gzip
import zipfile
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import copy, copy2
from datetime import date, datetime
import sys


sys.stderr = open(snakemake.log[0], "w")

file = str(snakemake.input)
name = os.path.splitext(file)[0]
shutil.copy(file, name + ".zip")
zipped_file = name + ".zip"

with zipfile.ZipFile(zipped_file, "r") as zip_ref:
    name = zipped_file.split("/")[-1]
    new_dir = str(snakemake.params.between) + "/" + name
    zip_ref.extractall(os.path.splitext(new_dir)[0] + "/")
directory = os.listdir(str(snakemake.params.between))
subdir = os.listdir(str(snakemake.params.between) + "/" + directory[0])
orig_dir = str(snakemake.params.between) + "/" + directory[0] + "/" + subdir[0]
new_dir = str(snakemake.params.between) + "/" + directory[0]
# print(new_dir)
dir_name = new_dir.split("/")[-1]
for f in os.listdir(orig_dir):
    path = orig_dir + "/" + f
    # print(path)
    unzipped_path = str(snakemake.params.between) + "/" + dir_name + "/" + f
    shutil.move(path, unzipped_path)
os.rmdir(orig_dir)

dirlist = os.listdir(str(snakemake.params.between))
subdir = []
i = 0
while i < len(dirlist):
    directory = str(snakemake.params.between) + "/" + dirlist[i] + "/"
    subdir.append(directory)
    i = i + 1

b = 0
while b < len(subdir):
    datadir = subdir[b] + "/" + "data/"
    if "human-count" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.human_count)
    b += 1
