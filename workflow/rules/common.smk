from snakemake.utils import validate
import pandas as pd
from datetime import datetime
from os import path, getcwd, listdir
from os.path import dirname
import gzip
import shutil
import os
import zipfile


configfile: "config/config.yaml"


def get_date():
    metadata = pd.read_csv(config["metadata"], header=0, delimiter=",")
    date = metadata["run_date"].iloc[1]
    return date


def get_samples():
    incoming_files = get_filenames()
    names = []
    for file in incoming_files:
        name = file.split("_")[0]
        names.append(name)
    return names


def get_filenames():
    # read files from sample_info
    data = pd.read_csv("config/pep/sample_info.txt")
    path_list1 = data["path1"].tolist()
    path_list2 = data["path2"].tolist()
    allpaths = path_list1 + path_list2
    incoming_files = []
    for file in allpaths:
        filename = file.split("/")[-1]
        if ".fastq.gz" in file:
            incoming_files.append(filename)
    return incoming_files


def get_data_dir():
    data = pd.read_csv("config/pep/sample_info.txt")
    path_list = data["path1"].tolist()
    try:
        path = path_list[0]
        dir = os.path.dirname(path)
        return dir
    except IndexError:
        print("There is no data directory jet known to the workflow.")


def get_file_dir(name):
    dir = get_data_dir + name
    return dir


def get_abundance(path):
    with open(path, "r") as f:
        abundance = f.read()
        return abundance
