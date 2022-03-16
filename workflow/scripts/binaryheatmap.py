import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
from pylab import savefig

# Reading information from a file, holding binary information whether an OTU is present or absent in a sample

# Reading the file, deleting unnecessary columns
df = pd.read_csv(str(snakemake.input), delimiter="\t", header=1, index_col="taxonomy")
del df["#OTU ID"]
# Creating a specific colour scheme
my_colors = [(0.2, 0.6, 0.3), (0.1, 0.4, 0.5)]
fig, ax = plt.subplots(figsize=(11, 9))
# Creating a heatmap with seaborn
heatmap = sb.heatmap(
    df, cmap=my_colors, linewidth=0.3, linecolor="Black", cbar=False, yticklabels=True
)
# Setting the properties for the legend
cbar = heatmap.figure.colorbar(heatmap.collections[0])
cbar.set_ticks([0.25, 0.75])
cbar.set_ticklabels(["absent", "present"])
heatmap.set_yticklabels(heatmap.get_ymajorticklabels(), fontsize=6)
heatmap.set(xlabel="Sample", ylabel="Taxonomic classification of bacteria")
# Setting the orientation of the plot
plt.subplots_adjust(bottom=0.2, left=0.7)
plt.show()

figure = heatmap.get_figure()
figure.savefig(str(snakemake.output), dpi=400)
