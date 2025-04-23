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
    if config["datatype"] == "SampleData[PairedEndSequencesWithQuality]":
        path_list1 = data["path1"].tolist()
        path_list2 = data["path2"].tolist()
        allpaths = path_list1 + path_list2
    elif config["datatype"] == "SampleData[SequencesWithQuality]":
        path_list1 = data["path1"].tolist()
        allpaths = path_list1
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
        direc = os.path.dirname(path)
        return direc
    except IndexError:
        print("There is no data directory jet known to the workflow.")


def get_file_dir(name):
    dir = get_data_dir + name
    return dir


def get_abundance(path):
    with open(path, "r") as f:
        abundance = f.read()
        return abundance


def get_for_testing():
    return bool(config["testing"])


def get_if_testing(string):
    return string if get_for_testing() else ""


def get_reads_for_kraken():
    incoming_files = get_filenames()
    names = []
    for file in incoming_files:
        name = file.split("_")[0]
        number = file.split("_")[1]
        kraken_name = "{name}_{number}".format(name=name, number=number)
        names.append(kraken_name)
        names = list(dict.fromkeys(names))
    return names


def get_metadata_columns():
    metadata = pd.read_csv(config["metadata"], header=0, delimiter=",")
    header = metadata.columns[metadata.isin(["numeric"]).any()].values.tolist()
    remove_list = config["remove-columns"]
    if len(remove_list) != 0:
        for item in remove_list:
            if item in header:
                try:
                    header.remove(item)
                except:
                    continue
    header = [item.lower() for item in header]
    return header


def get_metadata_categorical_columns():
    metadata = pd.read_csv(config["metadata"], header=0, delimiter=",")
    header = metadata.columns[metadata.isin(["categorical"]).any()].values.tolist()
    namelist = metadata.columns[metadata.isin(["categorical"]).any()].values.tolist()
    data_only = metadata.drop(labels=0, axis=0)
    y = 0
    while y < len(header):
        name = header[y]
        if (data_only[name] == data_only[name][1]).all():
            namelist.remove(name)
        if "barcode" in name:
            namelist.remove(name)
        y = y + 1
    remove_list = config["remove-columns"]
    if len(remove_list) != 0:
        for item in remove_list:
            if item in header:
                try:
                    namelist.remove(item)
                except:
                    continue
    namelist = [item.lower() for item in namelist]
    return namelist


def get_ancom_columns():
    return config["ancom"]["metadata-column"]


def get_metric(diversity):
    if diversity == "alpha":
        return config["diversity"]["alpha"]["diversity-metric"]
    if diversity == "beta":
        return config["diversity"]["beta"]["diversity-metric"]


def get_phylogenetic_metric(diversity):
    if diversity == "beta":
        return config["diversity"]["beta"]["phylogeny-metric"]
    if diversity == "alpha":
        return config["diversity"]["alpha"]["phylogeny-metric"]


def get_complete_beta_metric():
    complete = []
    complete.extend(config["diversity"]["beta"]["phylogeny-metric"])
    complete.extend(config["diversity"]["beta"]["diversity-metric"])
    return complete
