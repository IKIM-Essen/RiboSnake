from snakemake.utils import validate
import pandas as pd


##### load config and sample sheets #####

configfile: "config/config.yaml"


def get_samples():
    return list(pep.sample_table["sample_name"].values)