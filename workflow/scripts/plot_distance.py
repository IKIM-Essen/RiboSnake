import os
import shutil
import sys
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

sys.stderr = open(snakemake.log[0], "w")

directories = os.listdir(str(snakemake.input))

i = 0
while i < len(directories):
    output_name = os.path.split(str(snakemake.output))[1]
    name = output_name.split(".")[0]
    if name in directories[i]:
        df = pd.read_csv(
            str(snakemake.input)
            + "/"
            + directories[i]
            + "/data"
            + "/distance-matrix.tsv",
            sep="\t",
            header=0,
            index_col=0,
        )
        # Set the diagonal and lower triangle to NaN
        df.values[np.tril_indices_from(df)] = np.nan

        # Create a triangular matrix plot for the lower triangle
        fig = go.Figure(
            data=go.Heatmap(
                z=df.values,
                x=df.columns,
                y=df.index,
                colorscale="Viridis",  # Choose a colorscale
                zmin=0,  # Minimum value for color scale
                zmax=df.max().max(),  # Maximum value for color scale
                hoverongaps=False,  # Do not show hover information for NaN values
                colorbar=dict(
                    title="Distance",
                ),  # Add a title to the colorbar
            )
        )

        # Update the layout
        fig.update_layout(
            title="Beta diversity distance matrix",
            xaxis_title="Sample Names",
            yaxis_title="Sample Names",
        )

        # Save the plot as an HTML file
        fig.write_html(str(snakemake.output))
    i = i + 1
