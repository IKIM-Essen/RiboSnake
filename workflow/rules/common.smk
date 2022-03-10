from snakemake.utils import validate
import pandas as pd
from datetime import datetime
from os import path, getcwd, listdir
import gzip
import shutil
import os
import zipfile

##### load config and sample sheets #####


def get_date():
    metadata = pd.read_csv(config["metadata"], header = 0, delimiter = ",")
    date = metadata["run_date"].iloc[1]
    print(date)
    #date = datetime.today().strftime('%Y-%m-%d')
    return date

def get_samples():
    IN_PATH = get_data_dir()
    incoming_files = []
    for f in listdir(IN_PATH):
        if (
            path.isfile(path.join(IN_PATH, f))
            and ".fastq.gz" in f
        ):
            incoming_files.append(f)
    names = []
    for file in incoming_files:
        name  = file.split("_")[0]
        names.append(name)
    return names

def get_filenames():
    IN_PATH = get_data_dir()
    incoming_files = []
    for f in listdir(IN_PATH):
        if (
            path.isfile(path.join(IN_PATH, f))
            and ".fastq.gz" in f
        ):
            incoming_files.append(f)
    return incoming_files
    #for file in incoming_files:
    #    return file

def get_data_dir():
    dir = str(config["data"]) + get_date()
    #current_dir = dir + get_date()
    return dir

def get_file_dir(name):
    dir = get_data_dir+name
    return dir

def get_abundance(path):
    with open(path, "r") as f:
        abundance = f.read()
        return abundance

