import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from pylab import savefig


parser = argparse.ArgumentParser(
    description="Create interactive absolute taxa barplot"
)
parser.add_argument(
    "--input",
    required=True,
    help="Input taxonomy/feature count table"
)
parser.add_argument(
    "--output",
    required=True,
    help="Output HTML file"
)
parser.add_argument(
    "--samplename",
    required=True,
    help="Sample name metadata parameter"
)
parser.add_argument(
    "--metadata",
    required=True,
    help="Metadata TSV file"
)

args = parser.parse_args()
df = pd.read_csv(
    args.input,
    delimiter="\t",
    header=1,
    index_col="taxonomy"
)


if "#OTU ID" in df.columns:
    df.drop("#OTU ID", axis=1, inplace=True)
df_reduced = df.groupby(df.index).sum()
species = df_reduced.index.unique()
for i in df_reduced.index:
    if "Unassigned" not in i:
        df_reduced.index = (
            df_reduced.index
            .str.split("; s__")
            .str[0]
        )

df_reduced_second = df_reduced.groupby(
    df_reduced.index
).sum()
for i in list(df_reduced_second.index):
    if "Chloroplast" in i:
        df_reduced_second.drop(
            index=i,
            inplace=True
        )

df_reduced_second.replace(
    -np.inf,
    0,
    inplace=True
)
df_log10 = np.log10(
    df_reduced_second
)


metadata = pd.read_csv(
    args.metadata,
    delimiter="\t",
    header=0,
    index_col="sample-ID"
)

if "#q2:types" in metadata.index:
    metadata.drop(
        "#q2:types",
        axis=0,
        inplace=True
    )

dropdown_options = [
    {
        "label": col,
        "value": col
    }
    for col in metadata.columns
]
initial_x = metadata.columns[0]
x_labels = []
for sample_id in df_reduced_second.columns:

    if sample_id in metadata.index:
        x_labels.append(
            f"{sample_id}: "
            f"{metadata.at[sample_id, initial_x]}"
        )
    else:
        x_labels.append(
            str(sample_id)
        )
fig = go.Figure()

color_map = {}

for i, bacterium_name in enumerate(
    df_reduced_second.index
):
    color_map[bacterium_name] = (
        f"rgb("
        f"{i * 30 % 256}, "
        f"{i * 50 % 256}, "
        f"{i * 70 % 256}"
        f")"
    )


def add_traces(fig, x_labels):
    fig.data = []
    for bacterium_name, row in (
        df_reduced_second.iterrows()
    ):
        fig.add_trace(
            go.Bar(
                x=x_labels,
                y=row.values,
                name=bacterium_name,

                marker_color=[
                    color_map[bacterium_name]
                    for _ in x_labels
                ],
                marker=dict(
                    line=dict(width=0)
                ),
                hovertemplate=(
                    "<b>Bacterium:</b> "
                    + bacterium_name
                    + "<br>"
                    + "<b>Sample:</b> %{x}<br>"
                    + "<b>Abundance:</b> %{y:.2f}"
                    + "<br>"
                    + "<extra></extra>"
                )
            )
        )


add_traces(
    fig,
    x_labels
)
num_samples = len(
    df_reduced_second.columns
)
num_bacteria = len(
    df_reduced_second.index
)


min_width = 800

width = max(
    min_width,
    num_samples * 50
)

max_label_length = (
    max(len(label) for label in x_labels)
    if x_labels
    else 0
)

bottom_margin = max(
    150,
    int(max_label_length * 5)
)


max_bacteria_name_length = (
    max(
        len(name)
        for name in df_reduced_second.index
    )
    if len(df_reduced_second.index) > 0
    else 10
)
right_margin = max(
    250,
    int(max_bacteria_name_length * 6)
)


height = max(
    800,
    200 + num_bacteria * 8
)

buttons = []
for col in metadata.columns:
    button_labels = []
    for sample_id in df_reduced_second.columns:
        if sample_id in metadata.index:
            button_labels.append(
                f"{sample_id}: "
                f"{metadata.at[sample_id, col]}"
            )
        else:
            button_labels.append(
                str(sample_id)
            )
    buttons.append(
        {
            "args": [
                {
                    "x": [button_labels]
                }
            ],
            "label": col,
            "method": "restyle"
        }
    )
    
fig.update_layout(
    title=(
        "Barplot logarithmic "
        "absolute bacterial abundances"
    ),
    xaxis_title="Sample",
    yaxis_title=(
        "Logarithmic absolute "
        "bacterial abundance"
    ),
    barmode="stack",
    legend_title="Bacterial names",
    updatemenus=[
        {
            "buttons": buttons,
            "direction": "down",
            "showactive": True,
            "x": 1.02,
            "xanchor": "left",
            "y": 1.15,
            "yanchor": "top",
            "type": "dropdown"
        }
    ],
    legend=dict(
        orientation="v",
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.02,
        bgcolor="rgba(255, 255, 255, 0.8)",
        bordercolor="Black",
        borderwidth=1,
        font=dict(size=8),
        itemclick="toggleothers",
        traceorder="normal"
    ),
    xaxis=dict(
        tickangle=45,
        tickfont=dict(size=8)
    ),
    yaxis=dict(
        type="log",
        tickvals=[
            10 ** i
            for i in range(2, 15)
        ],
        ticktext=[
            f"10^{i}"
            for i in range(2, 15)
        ],
        tickmode="array",
        tickfont=dict(size=10)
    ),
    margin=dict(
        l=100,
        r=right_margin,
        t=50,
        b=bottom_margin
    ),
    width=width,
    height=height
)

fig.write_html(
    args.output
)
