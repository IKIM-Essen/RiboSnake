from pathlib import Path
import pandas as pd
from shutil import copy2
import gzip
import re
from collections import defaultdict
import datetime
import os

def validate_date(date_text):
        try:
            datetime.date.fromisoformat(date_text)
        except ValueError:
            raise ValueError("Incorrect data format, should be YYYY-MM-DD")

def print_errors(sample, errors):
    print(f"\nERROR: Sample '{sample}' failed validation:")
    for err in errors:
        print(f"   {err}")

def is_valid_fastq_gz(path):
    try:
        with gzip.open(path, "rt") as fh:
            for i in range(4):
                fh.readline()
        return True
    except Exception:
        return False

def count_reads_fastq_gz(path,MIN_READS):
    count = 0
    threshold = MIN_READS*4
    with gzip.open(path, "rt") as fh:
        for _ in fh:
            count += 1
            if count >= threshold:
                break
        if count < threshold:
            raise ValueError(f"{path} has fewer than {MIN_READS} reads")

def check_metadata_file(path):
    # Check q2 types line
    with open(path) as fh:
        lines = [l.strip() for l in fh if l.strip()]

    header = next((l for l in lines if l.startswith("sample_name")), None)
    q2_line = next((l for l in lines if l.startswith("#q2:types")), None)

    if q2_line is None:
        raise ValueError("Missing #q2:types line in metadata")

    # check every column namem has a q2 type
    q2_list = q2_line.strip().split(",")
    header_list = header.strip().split(",")
    if len(q2_list) != len(header_list):
        print(len(header_list),header_list)
        print(len(q2_list),q2_list)
        raise ValueError("Number of q2 types does not match metadata columns")

    required_columns = {"sample_name", "subject", "run_date"}
    missing = required_columns - set(metadata.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")
    
    if not 'run_date' in header_list:
        raise ValueError("metadata.txt is missing the mandatory run_date column")

    # Data format correct
    metadata_dates = list(set(metadata[metadata["run_date"] != 'categorical'].run_date))
    if len(metadata_dates) > 1:
        print(len(metadata_dates),metadata_dates)
        raise ValueError("Only one run_date acceptable")

    for date in metadata_dates:
        validate_date(date)

def validate_incoming(incoming,MIN_READS,FASTQ_PATTERN):
    parsed = []
    parse_errors = defaultdict(list)

    # check name formats correct
    for fq in incoming:
        m = FASTQ_PATTERN.match(fq.name)
        if not m:
            parse_errors[fq.name].append("Filename does not match required format")
            continue
        parsed.append((fq, m.groupdict()))
    # check every fastq determined in sample_name present in config["input"]
    fastq_samples = {p[1]["sample"] for p in parsed}

    # Missing FASTQs
    for sample in metadata_samples - fastq_samples:
        print_errors(sample, ["Sample defined in metadata but no FASTQ files found"])

    # Extra FASTQs
    for sample in fastq_samples - metadata_samples:
        print_errors(sample, ["FASTQ files present but sample not listed in metadata"])

    by_sample = defaultdict(list)
    for fq, info in parsed:
        by_sample[info["sample"]].append((fq, info))

    for sample, entries in by_sample.items():
        errors = []

        # check paired-end read presence
        reads = {e[1]["read"] for e in entries}
        if reads != {"1", "2"}:
            errors.append("Missing paired-end reads (R1/R2)")
        else: paired_end = True

        # check present gzip format
        for fq, _ in entries:
            if not is_valid_fastq_gz(fq):
                errors.append(f"{fq.name} is not a valid gzipped FASTQ")
                continue

            # has every read > 100 reads
            count_reads_fastq_gz(fq,MIN_READS)

        if errors:
            print_errors(sample, errors)

    # check metadata file run_date

    for fname, errs in parse_errors.items():
        print_errors(fname, errs)

sys.stderr = open(snakemake.log[0], "w")

# Creating a metadata sample-sheet, that holds additional information concerning the samples.
# Metadata information is needed to be able to create plots and for metadata-specific analysis.

config = snakemake.config

data_root = Path(os.path.abspath(config["data"]))
in_path = Path(os.path.abspath(config["input"]))
metadata = pd.read_csv(str(snakemake.input),comment='#', header=0, delimiter=",")
metadata.dropna(axis=1,inplace=True)
metadata_sample_list = metadata.sample_name.tolist()
metadata_samples = set(metadata.sample_name)

if len(metadata_sample_list) != len(metadata_samples):
    raise ValueError("Sample names not unique per sample")

MIN_READS = 100
FASTQ_PATTERN = re.compile(
    r"^(?P<sample>[A-Za-z0-9_.-]+)_(?P<num>\d+)_L(?P<lane>\d+)_R(?P<read>[12])_001\.fastq\.gz$"
)

# 0. Check metadata.to
check_metadata_file(str(snakemake.input))

# 1. collect FASTQ files
incoming = list(in_path.glob("*.fastq.gz"))

# 2. Validate FASTQ files
validate_incoming(incoming,MIN_READS,FASTQ_PATTERN)

# 3. read date and create date directory
date = metadata.loc[1, "run_date"]
data_path = data_root / date
data_path.mkdir(exist_ok=True)

# 4. copy files to date directory
existing = {f.name for f in data_path.iterdir() if f.is_file()}
to_copy = [f for f in incoming]

for f in to_copy:
    print(f"copying {f} to {data_path}")
    copy2(f, data_path)

# 5. sample names from filenames
sample_list = {f.name.split("_")[0] for f in incoming}

# 6. validate against metadata
sample_list = sample_list & metadata_samples

# 7. clean metadata
metadata = metadata.rename(columns={"sample_name": "sample-ID"}).fillna(0)
metadata.to_csv(snakemake.output.sample_tsv, sep="\t", index=False)

# 8. build path assignments
fastq_map = {s: {"R1": None, "R2": None} for s in sample_list}

for f in to_copy:
    sid = f.name.split("_")[0]
    new_path = str(f).replace("incoming",f"data/{date}")
    if sid in fastq_map:
        if "R1" in f.name:
            fastq_map[sid]["R1"] = str(new_path)
        if "R2" in f.name:
            fastq_map[sid]["R2"] = str(new_path)

sample_info = metadata.set_index("sample-ID")[["run_date"]].loc[list(sample_list)]

sample_info["path1"] = sample_info.index.map(lambda s: fastq_map[s]["R1"])
sample_info["path2"] = sample_info.index.map(lambda s: fastq_map[s]["R2"])

sample_info.to_csv(snakemake.output.sample_info)