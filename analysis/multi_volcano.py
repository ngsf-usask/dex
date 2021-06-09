# %%
import pandas as pd
import plotly.offline as pyo
import plotly.graph_objs as go
import numpy as np
import argparse

# LIMIT will never be able to select specific gene because limited by Scatter function. To employ Dash requires using a host server
# such as heroku


# %%
# Use command line functions
def get_args():
    """
    Purpose:
        Pull arguments from command line
    Return:
        arguments:
            .name : name of final graph
            .files : list of dataset name and file separated by space
                ex: dataset1 dataset1.csv
    """
    parser = argparse.ArgumentParser(description='Provide input for volcano plot')
    parser.add_argument('name', help='name of graph')
    parser.add_argument('files', help='file-of-file input, one on each line')

    arguments = parser.parse_args()
    return arguments

args = get_args()
graph_name = args.name


# Grab all the data and convert into dataframe
dfs = []
names = []

fof = args.files

with open(fof, 'r') as files:
    for file in files:
        # Take each line. Split, first item becomes the datasets name and the second one is stripped/cast and opened for a dataframe
        file = file.split(" ")
        names.append(file[0])
        file = str(file[1].strip())
        df = pd.read_csv(file, header=0)
        
        # Create x-axis values for log2 and then multiply by -1 if "DOWN" regulated
        df["log2(fold_change)"] = np.log2(df["fold_change"])
        df["log2(fold_change)"] = np.where(df["change_direction"] == "DOWN", df["log2(fold_change)"] * -1, df["log2(fold_change)"])

        # Create y-axis values for -log10 of the p value
        df["-log10(padj)"] = np.log10(df["padj"])*-1
        dfs.append(df)

# %%
# Create a list of every gene name in the DE files
gene_names = []
for df in dfs:
    for gene in df["gene_name"]:
        gene_names.append(gene)

# Remove duplicates, and sort
gene_names = list( dict.fromkeys(gene_names) )
gene_names.sort()


# %%
# create the graphs
data = []
count = 0

for df in dfs:
    trace = go.Scatter(
        x = df["log2(fold_change)"],
        y = df["-log10(padj)"],
        mode = 'markers',
        name= names[count],
        opacity = 0.8,
        hovertext= df["gene_name"],
        hovertemplate = "<b>Gene: %{hovertext}</b><br>" +
                        "<i>X</i>: %{x:.3f}<br>" +
                        "<i>Y</i>: %{y:.3f}",
        # hoverinfo='text',
        marker = dict(
            size = 12,
            line= dict(width=2,
                    color="DarkSlateGrey"),
            )
    )
    count += 1
    data.append(trace)

layout = go.Layout(
    title = graph_name, 
    xaxis = dict(title = "log2(fold_change)",
                # range = [-3.5,3.5],
                autorange=True, 
                ),
    yaxis = dict(title = "-log10(padj)",
                # range = [0,3.5],
                rangemode='tozero',
                ), 
    hovermode ='closest', 
)

fig = go.Figure(data=data, layout=layout)

# %%
if __name__ == "__main__":
    pyo.plot(fig, filename=graph_name + "_arg_volcano.html",
            )