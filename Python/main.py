import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from myplotting import *
import plotly
import plotly.io as pio


properties = ["massflow", "speed", "mu", "mubar", "dP", "dP_m",
              "demand", "power", "dPower", "TIn", "TOut", "dT", "delta", "nubar", "epsilon",
              "material_cost", "engineering_cost", "operating_cost", "resources_cost"]

# casedir = "C:/Users/cornelia.blanke/Documents/adv/Cases/Esta_NEW"
# case = "Esta_NEW"
# meteoyear = "2019"
casedir = "C:/Users/cornelia.blanke/Documents/adv/Cases/Projet"
case = "Projet1"
meteoyear = "2021"

df = {}     # empty dictionary for the dataframes


# read the files into dataframes
file = os.path.join(casedir, case + ".csv")
if os.path.exists(file):
    df["network"] = pd.read_csv(file, sep='\t', index_col=[0, 1])
    df["network"] = df["network"].rename(columns={'0.1': '0'})
# print(df["network"])
# fig = plot_tree(df["network"])


for prop in properties:
    file = os.path.join(casedir, case + "_results_" + prop + ".csv")
    if os.path.exists(file):
        df[prop] = pd.read_csv(file, sep='\t').set_index("Timestep")
    if prop in ["TIn", "TOut"]:
        df[prop] -= 273.15      # in Celsius


''' pictures for publication '''
publidir = "C:/Users/cornelia.blanke/OneDrive - HESSO/Cornelia/ADVENS/Publication"

pio.templates.default = 'simple_white'
#colors = px.colors.sequential.Greys
colors = ['rgb(115,115,115)', 'rgb(37,37,37)']
def save_fig(fig, figname):
    fig.update_layout(title=None,
                      autosize=False, width=470, height=300,
                      margin=dict(l=50, r=10, b=50, t=10, pad=0))
    try:
        for i, d in enumerate(fig.data):
            d.line.color = colors[i]
            d.line.width = i+1
    except:
        for i, d in enumerate(fig.data):
            d.marker.color = colors[i]
        fig.update_layout(width=700)
    fig.write_image(os.path.join(publidir, figname))

#fig = plot_over_time("massflow", "Total", df)
#fig = plot_over_time_sorted("massflow", "Total", df)
#fig = plot_over_time_cumulated("massflow", "Total", df)
#fig = plot_over_time_both("delta", "Central Plant", df, ascending=True)
#fig = plot_elements("mubar", "Branch", df)
#fig = plot_boxplot("TIn", "Segment", df)


G = compute_tree(df["network"])
fig = plot_tree(G)
fig.update_layout(title=None)
fig.data[0].line.color = 'Black'
fig.data[0].line.width = 2
fig.data[2].marker = {'color': 'White', 'size': 20, 'symbol': 'square', 'line': dict(color='Black', width=2)}
fig.data[3].marker = {'color': 'White', 'size': 12, 'symbol': 'circle', 'line': dict(color='Black', width=2)}
fig.data[4].marker = {'color': 'Black', 'size': 8, 'symbol': 'circle'}
#print(fig.data)
fig.write_image(os.path.join(publidir, case + "_tree.svg"))

df_meteo = read_meteo(os.path.join(casedir, case + "_meteo.csv"), meteoyear)
fig = plot_meteo(df_meteo, violin=False)
save_fig(fig, case + "_meteo.svg")

fig = plot_over_time_both("power", "Total", df, violin=False)
fig.update_layout(yaxis_title="Power [W]")
save_fig(fig, case + "_power.svg")

fig = plot_over_time("dPower", "Total", df, violin=False)
fig.update_layout(yaxis_title="Power Loss [W]")
save_fig(fig, case + "_dPower.svg")

fig = plot_over_time("TOut", "Total", df, violin=False)
fig.update_layout(yaxis_title="Temperature [°C]")
save_fig(fig, case + "_TOut.svg")

fig = plot_boxplot("TIn", "Substation", df)     # Plot Tmin?
fig.update_layout(yaxis_title="Temperature [°C]")
save_fig(fig, case + "_TIn_SST.svg")

if case == "Esta_NEW":
    fig = plot_over_time("TIn", "Substation: 23", df, violin=False)
    fig.update_layout(yaxis_title="Temperature [°C]")
    save_fig(fig, case + "_TIn_SST23.svg")

    fig.update_xaxes(range=[5000, 5500])
    save_fig(fig, case + "_TIn_SST23zoom.svg")

if case == "Projet1":
    fig = plot_over_time("TIn", "Substation: 11", df, violin=False)
    fig.update_layout(yaxis_title="Temperature [°C]")
    save_fig(fig, case + "_TIn_SST23.svg")

    fig = go.Figure([   # plot COP
        go.Scatter(x=df["demand"]["Substation: 11"].index,
                   y=df["demand"]["Substation: 11"] / (df["demand"]["Substation: 11"] - df["power"]["Substation: 11"]))
    ])
    fig.update_layout(xaxis=dict(domain=[0, 1], showline=True, mirror=True),
                      yaxis=dict(domain=[0, 1], showline=True, mirror='ticks'),
                      xaxis_title="Time [h]", yaxis_title="COP")
    save_fig(fig, case + "_COP_SST23.svg")


#plotly.offline.plot(fig)


