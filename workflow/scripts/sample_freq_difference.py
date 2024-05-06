import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

whuman = str(snakemake.params.visual_wh)
wohuman = str(snakemake.params.visual_woh)

whuman_df = pd.read_csv(whuman, sep=",", header=None, index_col=0)
wohuman_df = pd.read_csv(wohuman, sep=",", header=None, index_col=0)

whuman_df.columns = ["whuman"]
wohuman_df.columns = ["wohuman"]

combined = pd.concat([whuman_df, wohuman_df], axis=1)

for sample in combined.index:
    if combined.at[sample, "whuman"] == combined.at[sample, "wohuman"]:
        combined.drop([sample], inplace=True)
# combined.drop_duplicates(subset = ["whuman", "wohuman"], inplace = True)

for sample in combined.index:
    difference = combined.at[sample, "whuman"] - combined.at[sample, "wohuman"]
    combined.at[sample, "difference"] = difference
combined.index.name = "Sample"
combined.rename(
    columns={"whuman": "Reads with human", "wohuman": "Reads without human"},
    inplace=True,
)
combined.to_csv(str(snakemake.output))
