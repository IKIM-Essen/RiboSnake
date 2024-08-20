import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from pylab import savefig
import sys


sys.stderr = open(snakemake.log[0], "w")

df = pd.read_csv(str(snakemake.input), delimiter="\t", header=1, index_col="taxonomy")
df.drop("#OTU ID", axis=1, inplace=True)
df_reduced = df.groupby(df.index).sum()
species = df_reduced.index.unique()
for i in df_reduced.index:
    if "Unassigned" not in i:
        df_reduced.index = df_reduced.index.str.split("; s__").str[0]
df_reduced_second = df_reduced.groupby(df_reduced.index).sum()
for i in df_reduced_second.index:
    if "Chloroplast" in i:
        df_reduced_second.drop(index=i, inplace=True)
# Apply log10 transformation to the DataFrame
df_reduced_second.replace(-np.inf, 0, inplace=True)
df_log10 = np.log10(df_reduced_second)

metadata = pd.read_csv(
    str(snakemake.params.metadata), delimiter="\t", header=0, index_col="sample-ID"
)
metadata.drop("#q2:types", axis=0, inplace=True)

dropdown_options = [{"label": col, "value": col} for col in metadata.columns]
initial_x = metadata.columns[0]  # Set initial x-axis label

# Generate initial x-axis labels
x_labels = [
    f"{sample_id}: {metadata.at[sample_id, initial_x]}"
    for sample_id in df_reduced_second.columns
]

# Initialize the figure
fig = go.Figure()

color_map = {}
for i, bacterium_name in enumerate(df_reduced_second.index):
    color_map[bacterium_name] = f"rgb({i * 30 % 256}, {i * 50 % 256}, {i * 70 % 256})"


def add_traces(fig, x_labels):
    fig.data = []  # Clear existing traces
    for bacterium_name, row in df_reduced_second.iterrows():
        fig.add_trace(
            go.Bar(
                x=x_labels,  # x-axis labels from dropdown
                y=row.values,  # Values for the current bacterium
                name=bacterium_name,  # Bacterium name as legend label
                marker_color=[color_map[bacterium_name] for _ in x_labels],
                marker=dict(line=dict(width=0)),  # Remove the bar outline
            )
        )


# Add initial traces
add_traces(fig, x_labels)
"""
# Loop through each row in the DataFrame
for bacterium_name, row in df_reduced_second.iterrows():
    # Add a bar trace for each bacterium
    fig.add_trace(
        go.Bar(
            x=df_reduced_second.columns,  # Sample names on x-axis
            y=row.values,  # Values for the current bacterium
            name=bacterium_name,  # Bacterium name as legend label
            marker_color=[
                color_map[bacterium_name] for col in df_reduced_second.columns
            ],  # Color by sample name
            marker=dict(line=dict(width=0)),  # Remove the bar outline
        )
    )
"""

# Update layout for the plot
fig.update_layout(
    title="Barplot logarithmic absolute bacterial abundances",
    xaxis_title="Sample",
    yaxis_title="Logarithmic absolute bacterial abundance",
    barmode="stack",  # Stacked bar mode
    legend_title="Bacterial names",  # Legend title
    updatemenus=[
        {
            "buttons": [
                {
                    "args": [
                        {
                            "x": [
                                [
                                    f"{sample_id}: {metadata.at[sample_id, col]}"
                                    for sample_id in df_reduced_second.columns
                                ]
                            ]
                        }
                    ],
                    "label": col,
                    "method": "restyle",
                }
                for col in metadata.columns
            ],
            "direction": "down",
            "showactive": True,
            "x": 1.15,
            "xanchor": "left",
            "y": 1.2,
            "yanchor": "top",
            "type": "dropdown",
        }
    ],
    legend=dict(
        orientation="v",
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.01,
        bgcolor="rgba(255, 255, 255, 0.5)",  # Set legend background color with transparency
        bordercolor="Black",  # Set legend border color
        borderwidth=1,  # Set legend border width
        traceorder="normal",  # Set the order of legend items
    ),
    yaxis=dict(
        type="log",  # Set y-axis to logarithmic scale
        tickvals=[
            10 ** i for i in range(2, 15)
        ],  # Set tick values to 10^2, 10^3, and so on
        ticktext=[
            f"10^{i}" for i in range(2, 15)
        ],  # Set tick text to display 10^2, 10^3, and so on
        tickmode="array",
        tickfont=dict(size=10),  # Adjust tick font size
    ),
    margin=dict(
        l=100,  # Add left margin to accommodate tick text
        r=20,
        t=50,
        b=80,  # Add bottom margin to accommodate tick text
    ),
)
"""
# Update layout
fig.update_layout(
    title="Barplot logarithmic absolute bacterial abundances",
    xaxis_title="Sample",
    yaxis_title="Logarithmic absolute bacterial abundance",
    barmode="stack",  # Stacked bar mode
    legend_title="Bacterial names",
    updatemenus=[{
        "buttons": [
            {
                "args": [{"x": [metadata[col].values]}],
                "label": col,
                "method": "restyle"
            } for col in metadata.columns
        ],
        "direction": "down",
        "showactive": True,
        "x": 1.15,
        "xanchor": "left",
        "y": 1.2,
        "yanchor": "top",
        "type": "dropdown"
    }]
)
"""
# Save the figure as an HTML file
fig.write_html(str(snakemake.output))
