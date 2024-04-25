import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
import numpy as np
import plotly.graph_objects as go
from pylab import savefig
import sys


sys.stderr = open(snakemake.log[0], "w")

# Reading information from a file, holding binary information whether an OTU is present or absent in a sample

# Reading the file, deleting unnecessary columns
df = pd.read_csv(str(snakemake.input), delimiter="\t", header=1, index_col="taxonomy")
del df["#OTU ID"]

# Create color scheme
colorscale = [[0, "blue"], [1, "red"]]

# Create the heatmap figure
fig = go.Figure(
    data=go.Heatmap(z=df.values, x=df.columns, y=df.index, colorscale=colorscale)
)

# Update layout
fig.update_layout(
    title="Binary presence/absence heatmap",
    xaxis_title="Sample",
    yaxis_title="Taxon",
    height=800,  # Adjust height as needed
    width=1200,  # Adjust width as needed
)

# Save the figure as an HTML file
fig.write_html(str(snakemake.output))
# figure.savefig(str(snakemake.output), dpi=400)
