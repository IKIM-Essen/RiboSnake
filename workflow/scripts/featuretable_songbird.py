import pandas as pd
import os
import sys


sys.stderr = open(snakemake.log[0], "w")

differentials_df = pd.read_csv(
    str(snakemake.params), header=0, index_col="featureid", sep="\t"
)
taxa_df = pd.read_csv(str(snakemake.input.taxa), header=1, sep="\t")
taxa_small = taxa_df[["#OTU ID", "taxonomy"]].copy()
taxa_small_new = taxa_small.set_index(["#OTU ID"])
diff_taxa_df = pd.concat([differentials_df, taxa_small_new], axis=1)
diff_taxa_df_new = diff_taxa_df.set_index(["taxonomy"])
diff_taxa_df_new.to_csv(str(snakemake.output.diff), sep="\t", index=True)
taxa_small_new.rename_axis("featurue-id", inplace=True)
taxa_small_new.to_csv(str(snakemake.output.feature_meta), sep="\t", index=True)
