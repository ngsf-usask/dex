# %%
# Interactive python script created in vscode
# Uses plotly to build volcano plot overlaying multiple datasets

# TODO add ability to select specific gene

# %%
import pandas as pd
import plotly.offline as pyo
import plotly.graph_objs as go
import numpy as np

# %%
# Add name of data traces, file names, and final graph name
names = ["d4_vs_d0", "d1_vs_d0", "d4_vs_d1"]
files = ["experiment_d4_vs_d0_deseq2_DE_only.csv", "experiment_d1_vs_d0_deseq2_DE_only.csv",
            "experiment_d4_vs_d1_deseq2_DE_only.csv"]
graph_name = "all_days"

# Grab all the data and convert into dataframe
dfs = []
for file in files:
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


pyo.plot(fig, filename=graph_name + "_volcano.html",
        )
# %%
