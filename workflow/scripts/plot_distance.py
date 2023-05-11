import os
import shutil
import sys
import pandas as pd
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

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
        matrix = np.triu(df)
        fig, ax = plt.subplots(figsize=(20, 15))
        heatmap = sb.heatmap(
            df, mask=matrix, cbar_kws={"label": "Distance"}, annot=True, fmt=".1f"
        )
        heatmap.set(xlabel="Samplename", ylabel="Samplename")
        figure = heatmap.get_figure()
        figure.savefig(str(snakemake.output), dpi=400)
    i = i + 1
