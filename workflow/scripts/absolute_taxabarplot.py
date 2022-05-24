import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import savefig

df = pd.read_csv(str(snakemake.input), delimiter="\t", header=1, index_col="taxonomy")
print(df)
df.drop("#OTU ID", axis=1, inplace=True)
df_reduced = df.groupby(df.index).sum()
print(df_reduced)
species = df_reduced.index.unique()
for i in df_reduced.index:
    if "Unassigned" not in i:
        df_reduced.index = df_reduced.index.str.split("; s__").str[0]
df_reduced_second = df_reduced.groupby(df_reduced.index).sum()
for i in df_reduced_second.index:
    if "Chloroplast" in i:
        df_reduced_second.drop(index=i, inplace=True)
print(len(df_reduced_second.index))
df_trans = df_reduced_second.T
print(df_trans)
np.random.seed(100)
mycolors = np.random.choice(
    list(mpl.colors.XKCD_COLORS.keys()), len(species), replace=False
)
fig = df_trans.plot(kind="bar", stacked=True, color=mycolors, figsize=(16, 8))
handles, labels = fig.get_legend_handles_labels()
lgd = fig.legend(
    handles, labels, loc="center left", bbox_to_anchor=(1, 0.5), borderaxespad=0.0
)
plt.xlabel("Sample name")
plt.ylabel("Absolute species abundance")
plt.title("Taxa-bar-plot of absolute species abundance")
plt.savefig(str(snakemake.output), bbox_extra_artists=(lgd,), bbox_inches="tight")
