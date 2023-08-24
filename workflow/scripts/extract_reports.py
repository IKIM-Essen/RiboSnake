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
while b < len(subdir):
    datadir = subdir[b] + "data/"
    if "beta-rarefaction" in subdir[b]:
        svg = datadir + "heatmap.svg"
        shutil.copy(svg, snakemake.output.beta_svg)
        html = datadir
        shutil.copytree(html, snakemake.output.beta_html)
    elif "heatmap" in subdir[b] and "gneiss" not in subdir[b]:
        svg = datadir + "feature-table-heatmap.svg"
        shutil.copy(svg, snakemake.output.heatmap)
    elif "taxonomy" in subdir[b]:
        tsv = datadir + "metadata.tsv"
        shutil.copy(tsv, snakemake.output.taxonomy_tsv)
    elif "taxa-bar-plots" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.taxa_barplot)
    elif "heatmap_gneiss" in subdir[b]:
        svg = datadir + "heatmap.svg"
        shutil.copy(svg, snakemake.output.gneiss)
    if "alpha-rarefaction" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.alpha_html)
    if "paired-seqs" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.paired_seqs)
    if "fastq_stats" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.fastq_stats)
    if "demux-joined-filter-stats" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.demux_filter_stats)
    if "dada2-stats" in subdir[b]:
        html = datadir
        shutil.copytree(html, snakemake.output.dada2)
    b = b + 1
