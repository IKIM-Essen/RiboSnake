import pandas as pd

# Creating a csv file summarizing all important parameters
sys.stderr = open(snakemake.log[0], "w")


config = snakemake.config


Primer_error_rate = str(config["primertrimming"]["error_rate"])
Relative_abundance_filter = str(config["filtering"]["relative-abundance-filter"])
Min_seq_length = str(config["filtering"]["min-seq-length"])
Rarefaction_depth = str(config["rarefaction"]["sampling_depth"])
Taxa_heatmap_column = str(config["metadata-parameters"]["taxa-heatmap-column"])
Beta_metadata_column = str(config["metadata-parameters"]["beta-metadata-column"])
Gneiss_metadata_column = str(config["metadata-parameters"]["gneiss-metadata-column"])
Ancom_metadata_column = str(config["ancom"]["metadata-column"])

data = {
    "parameters":["Primer_error_rate", "Relative_abundance_filter", "Min_seq_length", "Rarefaction_depth", "Taxa_heatmap_column", "Beta_metadata_column", "Gneiss_metadata_column", "Ancom_metadata_column"],
    "values": [Primer_error_rate, Relative_abundance_filter, Min_seq_length, Rarefaction_depth, Taxa_heatmap_column, Beta_metadata_column, Gneiss_metadata_column, Ancom_metadata_column]
    }

df = pd.DataFrame(data = data)


print(df)
df.to_csv(str(snakemake.output), sep = ",", header = True, index = False)

