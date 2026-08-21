import argparse
import os
import pandas as pd
import shutil
import zipfile
import plotly.io as pio
from biom import load_table
import plotly.express as px

parser = argparse.ArgumentParser(
    description="Create rank abundance plots"
)
parser.add_argument(
    "--input",
    required=True,
    help="Input QZA file"
)
parser.add_argument(
    "--output-folder",
    required=True,
    help="Folder for individual plots"
)
parser.add_argument(
    "--output-file",
    required=True,
    help="Main HTML output file"
)
parser.add_argument(
    "--output-dir",
    required=True,
    help="Directory for temporary extraction"
)

args = parser.parse_args()


os.makedirs(
    args.output_folder,
    exist_ok=True
)
os.makedirs(
    args.output_dir,
    exist_ok=True
)
file = args.input
name = os.path.splitext(file)[0]
zip_file = name + ".zip"

shutil.copy(
    file,
    zip_file
)
with zipfile.ZipFile(
    zip_file,
    "r"
) as zip_ref:
    zip_name = os.path.basename(
        zip_file
    )
    extraction_base = os.path.join(
        args.output_dir,
        os.path.splitext(zip_name)[0]
    )
    os.makedirs(
        extraction_base,
        exist_ok=True
    )
    zip_ref.extractall(
        extraction_base
    )

directory_contents = os.listdir(
    extraction_base
)

if not directory_contents:
    raise RuntimeError(
        f"Keine Dateien im entpackten QZA-Verzeichnis: "
        f"{extraction_base}"
    )


data_dir = None
for item in directory_contents:

    candidate = os.path.join(
        extraction_base,
        item,
        "data"
    )

    if os.path.isdir(candidate):

        data_dir = candidate
        break


if data_dir is None:

    raise RuntimeError(
        "Konnte kein data/-Verzeichnis im "
        f"QZA finden: {extraction_base}"
    )



biom_file = os.path.join(
    data_dir,
    "feature-table.biom"
)

if not os.path.isfile(biom_file):

    raise FileNotFoundError(
        f"feature-table.biom wurde nicht gefunden: "
        f"{biom_file}"
    )


table = load_table(
    biom_file
)


df = pd.DataFrame(

    table.matrix_data.toarray(),

    columns=table.ids(
        axis="sample"
    ),

    index=table.ids(
        axis="observation"
    )
)



df = df.T

def create_interactive_plot(index):
    row = df.loc[index]
    sorted_abundance = (
        row.sort_values(
            ascending=False
        )
    )
    fig = px.bar(

        x=sorted_abundance.index,

        y=sorted_abundance.values,

        labels={
            "x": "Bacterial Names",
            "y": "Abundance"
        }
    )
    fig.update_layout(

        title=(
            f"Abundance of {index}"
        ),

        xaxis=dict(
            tickangle=-90,
            tickfont=dict(size=8)
        ),

        width=900,

        height=800,

        plot_bgcolor="white",

        paper_bgcolor="#f0f0f0"
    )

    return fig


#create all plots
all_figures = {

    index: create_interactive_plot(index)

    for index in df.index
}


for index, fig in all_figures.items():
    output_plot = os.path.join(

        args.output_folder,

        f"plot_{index}.html"
    )
    pio.write_html(
        fig,
        output_plot
    )

with open(
    args.output_file,
    "w"
) as f:

    f.write(
        "<html>\n"
    )

    f.write(
        "<head>\n"
    )

    f.write(
        "</head>\n"
    )

    f.write(
        "<body>\n"
    )

    f.write(
        '<select id="plot_selector" '
        'name="plot_selector" '
        'onchange="update_plot()">\n'
    )

    for index in df.index:
        f.write(
            f'<option value="{index}">'
            f'{index}'
            f'</option>\n'
        )
    f.write(
        "</select>\n"
    )
    f.write(
        '<div id="plot_area">\n'
    )
    f.write(
        "</div>\n"
    )
    f.write(
        "<script>\n"
    )
    f.write(
        "function update_plot() {\n"
    )
    f.write(
        '  var selector = '
        'document.getElementById("plot_selector");\n'
    )
    f.write(
        "  var index = selector.value;\n"
    )
    f.write(
        '  var plot_div = '
        'document.getElementById("plot_area");\n'
    )
    f.write(
        "  plot_div.innerHTML = "
        "'<iframe src=\"plot_' + index + "
        "'.html\" width=\"1000\" "
        "height=\"800\" "
        "frameborder=\"0\"></iframe>';\n"
    )
    f.write(
        "}\n"
    )
    f.write(
        "</script>\n"
    )
    f.write(
        "</body>\n"
    )
    f.write(
        "</html>\n"
    )