import qiime2
import pandas as pd
import gzip
import shutil
import os
import zipfile
import sys


sys.stderr = open(snakemake.log[0], "w")

# Load the Qiime2 feature table artifact
feature_table = qiime2.Artifact.load(str(snakemake.input))

# Convert the feature table to a Pandas DataFrame
table_df = feature_table.view(pd.DataFrame)

# Function to calculate relative abundance for each sample
def calculate_relative_abundance(df):
    return df.div(df.sum(axis=1), axis=0)

# Calculate relative abundances for each sample
relative_abundance_df = calculate_relative_abundance(table_df)
print(relative_abundance_df)

# Set the minimum threshold for relative abundance in each sample
min_abundance = float(str(snakemake.params))  # Adjust this as needed

# Filter out features in each sample that are below the threshold
filtered_df = relative_abundance_df.applymap(lambda x: x if x >= min_abundance else 0)

# Filter out rows (features) that are all zeros across samples
filtered_df = filtered_df.loc[~(filtered_df == 0).all(axis=1)]

# Convert back to the original count format (for Qiime2 compatibility)
filtered_count_df = (filtered_df.mul(table_df.sum(axis=1), axis=0)).astype(int)

# Convert the filtered DataFrame back to a Qiime2 Artifact
filtered_table_artifact = qiime2.Artifact.import_data('FeatureTable[Frequency]', filtered_count_df)

# Save the filtered feature table as a Qiime2 artifact
filtered_table_artifact.save(str(snakemake.output))

print("Filtered feature table saved.")
