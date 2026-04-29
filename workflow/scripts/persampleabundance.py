import argparse
import qiime2
import pandas as pd
import sys


def calculate_relative_abundance(df):
    return df.div(df.sum(axis=1), axis=0)


def filter_feature_table(input_table, output_table, min_abundance, log_path=None):
    if log_path is not None:
        sys.stderr = open(log_path, "w")

    feature_table = qiime2.Artifact.load(str(input_table))
    table_df = feature_table.view(pd.DataFrame)

    relative_abundance_df = calculate_relative_abundance(table_df)
    min_abundance = float(min_abundance)

    filtered_df = relative_abundance_df.applymap(lambda x: x if x >= min_abundance else 0)
    filtered_df = filtered_df.loc[~(filtered_df == 0).all(axis=1)]
    filtered_count_df = (filtered_df.mul(table_df.sum(axis=1), axis=0)).astype(int)

    filtered_table_artifact = qiime2.Artifact.import_data("FeatureTable[Frequency]", filtered_count_df)
    filtered_table_artifact.save(str(output_table))

    print("Filtered feature table saved.")


def main():
    parser = argparse.ArgumentParser(description="Filter a Qiime2 feature table by per-sample relative abundance.")
    parser.add_argument("input_table", help="Input feature table artifact path")
    parser.add_argument("min_abundance", help="Minimum relative abundance threshold")
    parser.add_argument("output_table", help="Output filtered feature table artifact path")
    parser.add_argument("--log", help="Optional log file path", default=None)
    args = parser.parse_args()

    filter_feature_table(args.input_table, args.output_table, args.min_abundance, args.log)


if __name__ == "__main__":
    main()
else:
    filter_feature_table(
        snakemake.input.table,
        snakemake.output.table,
        snakemake.params.min_abundance,
        snakemake.log[0],
    )
