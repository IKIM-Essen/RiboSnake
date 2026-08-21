import pandas as pd
import argparse


parser = argparse.ArgumentParser(
    description="Calculate differences between read frequencies with and without human reads."
)

parser.add_argument(
    "--whuman",
    required=True,
    help="CSV file containing frequencies with human reads"
)

parser.add_argument(
    "--wohuman",
    required=True,
    help="CSV file containing frequencies without human reads"
)

parser.add_argument(
    "--output",
    required=True,
    help="Output CSV file"
)

args = parser.parse_args()



# Read input tables

whuman_df = pd.read_csv(
    args.whuman,
    sep=",",
    header=None,
    index_col=0
)

wohuman_df = pd.read_csv(
    args.wohuman,
    sep=",",
    header=None,
    index_col=0
)



# Rename columns

whuman_df.columns = ["whuman"]
wohuman_df.columns = ["wohuman"]


# Combine tables

combined = pd.concat(
    [whuman_df, wohuman_df],
    axis=1
)


# Remove samples where read counts are identical

combined = combined[
    combined["whuman"] != combined["wohuman"]
].copy()


combined["difference"] = (
    combined["whuman"] - combined["wohuman"]
)


combined.index.name = "Sample"

combined.rename(
    columns={
        "whuman": "Reads with human",
        "wohuman": "Reads without human"
    },
    inplace=True
)


combined.to_csv(args.output)