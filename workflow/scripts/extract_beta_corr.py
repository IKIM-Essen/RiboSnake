import os
import pandas as pd
import shutil
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import copy, copy2
from datetime import date, datetime
import sys


sys.stderr = open(snakemake.log[0], "w")
# Function to extract results from the qiime2 artifacts, already saved in a zip file

"""print("results/2022-02-24/visual/unzipped/")
os.walk("results/2022-02-24/visual/unzipped/")
for x in os.walk("results/2022-02-24/visual/unzipped/"):
    print(x[0])"""
# Reading the zip-files from the directory
dirlist = os.listdir(str(snakemake.input))
# Getting the subdirectories of the files
subdir = []
i = 0
while i < len(dirlist):
    directory = str(snakemake.input) + "/" + dirlist[i] + "/"
    subdir.append(directory)
    i = i + 1
"""z = 0
while z < len(subdir):
    if path.isdir(subdir[z]):
        dir = os.listdir(subdir[z])
        print(dir)
        subdir[z] = subdir[z] + dir[0] + "/"
        print(subdir[z])
        z = 1 + z
    else:
        z = 1 + z"""
# Iterating through the different subdirectories and copying the important information into the output directory
b = 0
path_name = str(snakemake.output)
# print(path_name)
plot_name = path_name.split("/")[-1]
while b < len(subdir):
    datadir = subdir[b] + "data/"
    if plot_name == subdir[b].split("/")[-2]:
        html = datadir
        shutil.copytree(html, str(snakemake.output))
    b = b + 1
