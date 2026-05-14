import pandas as pd
import os
import matplotlib.pyplot as plt
import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def process_file(tsv_file, output_dir, suffix):
    logging.info(f"Processing {suffix} file: {tsv_file}")
    # Read the TSV file, skipping the comment line
    df = pd.read_csv(tsv_file, sep='\t', header = 1)
    logging.info(f"Loaded dataframe with shape: {df.shape}")
    print(df)
    
    # The columns are: #OTU ID, sample1, sample2, ..., taxonomy
    # Last column is taxonomy
    samples = df.columns[1:-1]  # Skip #OTU ID and taxonomy
    logging.info(f"Identified {len(samples)} samples: {list(samples)}")
    
    for sample in samples:
        logging.info(f"Processing sample: {sample}")
        # Create dataframe with OTU ID, abundance, taxonomy
        sample_df = df[['#OTU ID', sample, 'taxonomy']].copy()
        # Filter out zero abundances
        sample_df = sample_df[sample_df[sample] > 0]
        # Sort by abundance descending
        sample_df = sample_df.sort_values(sample, ascending=False)
        logging.info(f"Sample {sample} has {len(sample_df)} non-zero taxa")
        
        # Save to TSV
        output_file = os.path.join(output_dir, f"{sample}_{suffix}.tsv")
        sample_df.to_csv(output_file, sep='\t', index=False)
        logging.info(f"Saved TSV: {output_file}")

# Main execution
output_dir = snakemake.output[0]
os.makedirs(output_dir, exist_ok=True)
logging.info(f"Output directory: {output_dir}")

process_file(snakemake.input.abs, output_dir, 'absolute')
process_file(snakemake.input.rel, output_dir, 'relative')

logging.info("Script completed successfully")