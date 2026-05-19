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
                hovertemplate="<b>Bacterium:</b> " 
                + bacterium_name + "<br>" 
                + "<b>Sample:</b> %{x}<br>" 
                + "<b>Abundance (log10):</b> %{y:.2f}<br>" 
                + "<extra></extra>",  # Remove secondary box
            )
        )


# Add initial traces
add_traces(fig, x_labels)

# Calculate dynamic margins and sizes based on data
num_samples = len(df_reduced_second.columns)
num_bacteria = len(df_reduced_second.index)

# Ensure minimum width for single or few samples
min_width = 800
width = max(min_width, num_samples * 50)

# Calculate bottom margin based on label length
max_label_length = max(len(label) for label in x_labels) if x_labels else 0
bottom_margin = max(150, int(max_label_length * 5))

# Calculate right margin for legend based on number of bacteria and label length
max_bacteria_name_length = (
    max(len(name) for name in df_reduced_second.index) 
    if len(df_reduced_second.index) > 0 
    else 10
)
right_margin = max(250, int(max_bacteria_name_length * 6))

# Calculate height based on number of bacteria
height = max(800, 200 + num_bacteria * 8)

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
            "x": 1.02,
            "xanchor": "left",
            "y": 1.15,
            "yanchor": "top",
            "type": "dropdown",
        }
    ],
    legend=dict(
        orientation="v",
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.02,
        bgcolor="rgba(255, 255, 255, 0.8)",  # Set legend background color with transparency
        bordercolor="Black",  # Set legend border color
        borderwidth=1,  # Set legend border width
        font=dict(size=8),
        itemclick="toggleothers",
        traceorder="normal",  # Set the order of legend items
    ),
    xaxis=dict(
        tickangle=45,  # Rotate x-axis tick labels by 45 degrees
        tickfont=dict(size=8),  # Adjust font size for x-axis labels
    ),
    yaxis=dict(
        type="log",  # Set y-axis to logarithmic scale
        tickvals=[
            10**i for i in range(2, 15)
        ],  # Set tick values to 10^2, 10^3, and so on
        ticktext=[
            f"10^{i}" for i in range(2, 15)
        ],  # Set tick text to display 10^2, 10^3, and so on
        tickmode="array",
        tickfont=dict(size=10),  # Adjust tick font size
    ),
    margin=dict(
        l=100,  # Add left margin to accommodate tick text
        r=right_margin,
        t=50,
        b=bottom_margin,  # Dynamic bottom margin to accommodate tick text
    ),
    width=width,  # Dynamic figure width based on sample count
    height=height,  # Dynamic figure height based on bacteria count
)

# Save the figure as an HTML file
fig.write_html(str(snakemake.output))
