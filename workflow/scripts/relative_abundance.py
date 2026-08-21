import argparse
import pandas as pd
import shutil
import os
import zipfile

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True)
parser.add_argument(
    "--output-abundance",
    required=True
)
parser.add_argument(
    "--feature-table",
    required=True
)
parser.add_argument(
    "--relative-abundance",
    required=True,
    type=float
)
args = parser.parse_args()
file = args.input

name = os.path.splitext(file)[0]
# Create temporary ZIP copy of the QZV
shutil.copy(file, name + ".zip")
filename = name + ".zip"

with zipfile.ZipFile(filename, "r") as zip_ref:
    name = filename.split("/")[-1]
    dir_name = os.path.dirname(file)
    new_dir = os.path.join(dir_name, name)
    # Remove previous extraction if it exists
    if os.path.isdir(new_dir):
        shutil.rmtree(new_dir)

    zip_ref.extractall(
        os.path.splitext(new_dir)[0] + "/"
    )

# Remove .zip from directory name
name = name.split(".")[0]
directory = os.path.join(
    os.path.dirname(file),
    name
)

# Move files from extracted subdirectories
for item in os.listdir(directory):
    orig_dir = os.path.join(
        directory,
        item
    )
    # Ignore files such as VERSION, metadata.yaml, etc.
    if not os.path.isdir(orig_dir):
        continue
    for f in os.listdir(orig_dir):
        path = os.path.join(
            orig_dir,
            f
        )
        destination = os.path.join(
            directory,
            f
        )
        # Only move files if the destination does not
        # already exist
        if not os.path.exists(destination):
            shutil.move(
                path,
                destination
            )


# Read sample-frequency-detail.csv
datadir = args.feature_table + "/"
csv = os.path.join(
    datadir,
    "sample-frequency-detail.csv"
)

if not os.path.exists(csv):
    raise FileNotFoundError(
        f"sample-frequency-detail.csv not found: {csv}"
    )

frequency = pd.read_csv(
    csv,
    header=None,
    delimiter=","
)
frequency.columns = [
    "Sample",
    "Abundance"
]

# Calculate abundance threshold
abundance = args.relative_abundance

column_sums = frequency.sum()


median_of_sums = column_sums.median()

endnumber = median_of_sums * abundance
endnumber = int(endnumber)
# Minimum threshold = 1
if endnumber == 0:
    endnumber += 1


# Write abundance threshold
with open(
    args.output_abundance,
    "w"
) as f:

    f.write(
        str(endnumber)
    )