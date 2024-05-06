import os
import pandas as pd
import shutil
import zipfile
from os.path import isfile, join
from os import path, getcwd, listdir
from os import mkdir
from shutil import copy, copy2
from datetime import date, datetime
import sys
import plotly.io as pio
from biom import load_table
import plotly.express as px

sys.stderr = open(snakemake.log[0], "w")

file = str(snakemake.input)
name = os.path.splitext(file)[0]
shutil.copy(file, name + ".zip")
zipped_file = name + ".zip"

with zipfile.ZipFile(zipped_file, "r") as zip_ref:
    name = zipped_file.split("/")[-1]
    new_dir = str(snakemake.params) + name
    zip_ref.extractall(os.path.splitext(new_dir)[0] + "/")

print(new_dir)
direc = new_dir.split(".")[0]
directory_contents = os.listdir(direc)
print(directory_contents)

datadir = new_dir.split(".")[0] + "/" + directory_contents[0] + "/data/"
print(datadir)
# Load the biom table
table = load_table(datadir + "feature-table.biom")

print(table.matrix_data.toarray())
print(table.ids(axis="observation"))
print(table.ids(axis="sample"))
# Convert biom table to pandas DataFrame
df = pd.DataFrame(
    table.matrix_data.toarray(),
    columns=table.ids(axis="sample"),
    index=table.ids(axis="observation"),
)
df = df.T
print(df)

# Function to create an interactive plot for each bacterial species
def create_interactive_plot(index):
    # Get the row corresponding to the index
    row = df.loc[index]

    # Sort the abundance values in descending order to determine plot order
    sorted_abundance = row.sort_values(ascending=False)

    # Create Plotly figure using Plotly Express
    fig = px.bar(
        x=sorted_abundance.index,
        y=sorted_abundance.values,
        labels={"x": "Bacterial Names", "y": "Abundance"},
    )

    # Update layout
    fig.update_layout(
        title=f"Abundance of {index}",
        xaxis=dict(tickangle=-90, tickfont=dict(size=8)),  # Set font size to 10
        width=900,  # Adjust width as needed
        height=800,  # Adjust height as needed
        plot_bgcolor="white",  # Set plot background color to white
        paper_bgcolor="#f0f0f0",  # Set paper background color to light gray
    )

    return fig


# Create a dictionary to store all the plotly figures
all_figures = {index: create_interactive_plot(index) for index in df.index}

# Save all the figures to HTML files
for index, fig in all_figures.items():
    pio.write_html(fig, str(snakemake.output.folder) + f"/plot_{index}.html")

# Create the HTML file with the dropdown selector and plot
with open(str(snakemake.output.file), "w") as f:
    f.write("<html>\n<head>\n</head>\n<body>\n")
    f.write(
        '<select id="plot_selector" name="plot_selector" onchange="update_plot()">\n'
    )
    for index in df.index:
        f.write(f'<option value="{index}">{index}</option>\n')
    f.write("</select>\n")
    f.write('<div id="plot_area">\n')
    f.write("</div>\n")
    f.write("<script>\n")
    f.write("function update_plot() {\n")
    f.write('  var selector = document.getElementById("plot_selector");\n')
    f.write("  var index = selector.value;\n")
    f.write('  var plot_div = document.getElementById("plot_area");\n')
    f.write(
        '  plot_div.innerHTML = \'<iframe src="plot_\' + index + \'.html" width="1000" height="800" frameborder="0"></iframe>\';\n'
    )
    f.write("}\n")
    f.write("</script>\n")
    f.write("</body>\n</html>\n")
