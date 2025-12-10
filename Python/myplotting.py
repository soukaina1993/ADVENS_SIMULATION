# import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
import plotly.express as px
import plotly.graph_objs as go
from plotly.subplots import make_subplots


def celsius(temperature):
    if temperature > 200:
        temperature -= 273.15
    return temperature


def read_meteo(meteofile, meteoyear, sep='\t'):
    if not isinstance(meteoyear, list):
        meteoyear = [meteoyear]
    df = pd.read_csv(meteofile, sep=sep, index_col='Timestep', usecols=['Timestep'] + meteoyear)
    df = df.applymap(celsius)
    return df


def plot_meteo(df, violin=True, colsize=0.9, spacing=0.01):
    x = list(df.index)
    y = df.T.values
    name = list(df.columns)
    fig = plot_widget(x, y, name=name, violin=violin, colsize=colsize, spacing=spacing)
    fig.update_layout(xaxis=dict(title_text='Time [h]'),
                      yaxis=dict(title_text='Temperature [°C]'))
    return fig


def plot_over_time(property, element, df_list, name=None, violin=True, colsize=0.9, spacing=0.01):
    if not isinstance(df_list, list):
        df_list = [df_list]
    x = df_list[0][property][element].index
    y = []
    for i in range(len(df_list)):
        if (df_list[i][property][element].index == x).all():
            y.append(df_list[i][property][element].values)
    title = element + " - " + property
    #fig = px.line(x=x, y=y, title=title)
    fig = plot_widget(x, y, name=name, violin=violin, colsize=colsize, spacing=spacing)
    n_color = len(fig.layout.template.layout.colorway)
    fig.update_traces(selector=dict(type='scatter'),
                      hovertemplate="Time = %{x}<br>" + property + " = %{y}" + '<extra></extra>')
    for i, d in enumerate(fig.data):
        if violin:
            idx = int(i/2)
        else:
            idx = i
        d['line']['color'] = fig.layout.template.layout.colorway[(3*idx) % n_color]
    fig.update_layout(title=dict(text=title, y=0.9),
                      xaxis_title="Time [h]", yaxis_title=None)
    return fig


def plot_over_time_sorted(property, element, df_list, name=None, violin=True, colsize=0.9, spacing=0.01):
    if not isinstance(df_list, list):
        df_list = [df_list]
    x = df_list[0][property][element].index
    y = []
    for i in range(len(df_list)):
        if (df_list[i][property][element].index == x).all():
            y.append(df_list[i][property][element].sort_values(ascending=False, ignore_index=True).values)
    title = element + " - " + property + " (sorted)"
    #fig = px.line(x=x, y=y, title=title)
    fig = plot_widget(x, y, name=name, violin=violin, colsize=colsize, spacing=spacing)
    n_color = len(fig.layout.template.layout.colorway)
    fig.update_traces(selector=dict(type='scatter'),
                      hovertemplate="Time = %{x}<br>" + property + " = %{y}" + '<extra></extra>')
    for i, d in enumerate(fig.data):
        if violin:
            idx = int(i/2)
        else:
            idx = i
        d['line']['color'] = fig.layout.template.layout.colorway[(3*idx+1) % n_color]
    fig.update_layout(title=dict(text=title, y=0.9), xaxis_title="Time [h]", yaxis_title=None)
    return fig


def plot_over_time_cumulated(property, element, df_list, name=None, violin=True, colsize=0.9, spacing=0.01):
    if not isinstance(df_list, list):
        df_list = [df_list]
    x = df_list[0][property][element].index
    y = []
    for i in range(len(df_list)):
        if (df_list[i][property][element].index == x).all():
            y.append(df_list[i][property][element].cumsum().values)
    title = element + " - " + property + " (cumulated)"
    #fig = px.line(x=x, y=y, title=title)
    fig = plot_widget(x, y, name=name, violin=violin, colsize=colsize, spacing=spacing)
    n_color = len(fig.layout.template.layout.colorway)
    fig.update_traces(selector=dict(type='scatter'),
                      hovertemplate="Time = %{x}<br>" + property + " = %{y}" + '<extra></extra>')
    for i, d in enumerate(fig.data):
        if violin:
            idx = int(i/2)
        else:
            idx = i
        d['line']['color'] = fig.layout.template.layout.colorway[(3*idx+2) % n_color]
    fig.update_layout(title=dict(text=title, y=0.9), xaxis_title="Time [h]", yaxis_title=None)
    return fig


def plot_over_time_both(property, element, df_list, name=None, ascending=False, violin=True, colsize=0.9, spacing=0.01):
    if not isinstance(df_list, list):
        df_list = [df_list]
    x = df_list[0][property][element].index
    y = []
    y_cum = []
    for i in range(len(df_list)):
        if (df_list[i][property][element].index == x).all():
            y.append(df_list[i][property][element].values)
            y_cum.append(df_list[i][property][element].sort_values(ascending=ascending, ignore_index=True).values)
    title = element + " - " + property
    # fig = px.line(x=x, y=[y, y_cum], title=title,
    #               color_discrete_sequence=['blue', 'red'])
    # newname = {'wide_variable_0': 'time-dependent', 'wide_variable_1': 'sorted'}
    # fig.for_each_trace(lambda t: t.update(name=newname[t.name]))
    if name is None:
        name_t = 'time-dependent'
        name_s = ['sorted'] * len(y_cum)
    else:
        name_t = [n + ' (time-dep)' for n in name]
        name_s = [n + ' (sorted)' for n in name]
    fig = plot_widget(x, y, name=name_t, violin=violin, colsize=colsize, spacing=spacing)
    n_color = len(fig.layout.template.layout.colorway)
    for i, d in enumerate(fig.data):
        if violin:
            idx = int(i/2)
        else:
            idx = i
        d['line']['color'] = fig.layout.template.layout.colorway[(3*idx) % n_color]
    # add sorted data
    for i in range(len(y_cum)):
        fig.add_trace(go.Scatter(x=x, y=y_cum[i], mode='lines',
                                 line=dict(color=fig.layout.template.layout.colorway[(3*i+1) % n_color]),
                                 name=name_s[i]))
    fig.update_traces(selector=dict(type='scatter'),
                      hovertemplate="Time = %{x}<br>" + property + " = %{y}" + '<extra></extra>')
    fig.update_layout(title=dict(text=title, y=0.9), xaxis_title="Time [h]", yaxis_title=None, legend_title=None)
    return fig


''' matplotlib '''
# def plot_elements(property, elementtype, df):
#     fig = plt.figure()
#     x = df[property].columns[df[property].columns.str.contains(elementtype)]
#     y = df[property][x]
#     ticks = range(len(x))
#     labels = x.str.replace(elementtype + ': ', '')
#     plt.stem(x, y.max(), 'tab:blue', basefmt=' ')
#     plt.stem(x, y.mean(), 'tab:orange', basefmt=' ')
#     plt.stem(x, y.min(), 'tab:green', basefmt=' ')
#     plt.xticks(ticks=ticks, labels=labels, rotation=90)
#     plt.xlabel(elementtype)
#     plt.title(elementtype + " - " + property)
#     return fig


def plot_elements(prop, elementtype, df, top=None, method='mean'):
    # top > 0 means "display n=top best elements regarding 'method'"
    # top < 0 means "display n=top worst elements regarding 'method'"
    x = df[prop].columns[df[prop].columns.str.contains(elementtype)]
    y = df[prop][x]
    if top is not None:
        df_agg = getattr(y, method)()
        if top > 0:
            x = df_agg.nlargest(n=top).index
        else:
            x = df_agg.nsmallest(n=-top).index
        y = df[prop][x]

    title = elementtype + " - " + prop
    if y.shape[0] > 1:
        data = [go.Scatter(x=x,  y=y.max(), mode='markers', marker=dict(color='blue'), name='max'),
                go.Scatter(x=x, y=y.mean(), mode='markers', marker=dict(color='red'), name='mean'),
                go.Scatter(x=x,  y=y.min(), mode='markers', marker=dict(color='darkgreen'), name='min')
                ]
        layout = go.Layout(
            shapes=[dict(type='line', x0=i, y0=0, x1=i, y1=y.max()[i],      # xref='x', yref='y',
                         line=dict(color='blue', width=1)) for i in range(len(x))] +
                   [dict(type='line', x0=i, y0=0, x1=i, y1=y.mean()[i],     # xref='x', yref='y',
                         line=dict(color='red', width=1)) for i in range(len(x))] +
                   [dict(type='line', x0=i, y0=0, x1=i, y1=y.min()[i],      # xref='x', yref='y',
                         line=dict(color='darkgreen', width=1)) for i in range(len(x))],
            title=title,
            xaxis=dict(tickmode='array', tickangle=90, tickvals=x,
                       ticktext=x.str.replace(elementtype + ': ', '')),
            legend=dict(yanchor="top", y=1, xanchor="right", x=1,
                        bgcolor='rgba(255,255,255,0.5)')
        )
    else:
        data = [go.Scatter(x=x, y=y.iloc[0, :], mode='markers', marker=dict(color='blue'))]
        layout = go.Layout(
            shapes=[dict(type='line', x0=i, y0=0, x1=i, y1=y.iloc[0, i],  # xref='x', yref='y',
                         line=dict(color='blue', width=1)) for i in range(len(x))],
            title=title,
            xaxis=dict(tickmode='array', tickangle=90, tickvals=x,
                       ticktext=x.str.replace(elementtype + ': ', '')),
        )
    fig = go.Figure(data, layout)
    fig.update_traces(hovertemplate=elementtype + ": %{x}<br>" + prop + " = %{y}" + '<extra></extra>')
    fig.update_layout(xaxis_title=elementtype, yaxis_title=None)
    return fig


''' seaborn '''
# def plot_boxplot(property, elementtype, df):
#     fig = plt.figure()
#     x = df[property].columns[df[property].columns.str.contains(elementtype)]
#     y = df[property][x]
#     ticks = range(len(x))
#     labels = x.str.replace(elementtype + ': ', '')
#     sns.boxplot(y, color='tab:cyan', width=0.5)
#     plt.xticks(ticks=ticks, labels=labels, rotation=90)
#     plt.xlabel(elementtype)
#     plt.title(elementtype + " - " + property)
#     return fig


def plot_boxplot(prop, elementtype, df, top=None, method='mean'):
    x = df[prop].columns[df[prop].columns.str.contains(elementtype)]
    y = df[prop][x]
    if top is not None:
        df_agg = getattr(y, method)()
        if top > 0:
            x = df_agg.nlargest(n=top).index
        else:
            x = df_agg.nsmallest(n=-top).index
        y = df[prop][x]

    title = elementtype + " - " + prop
    fig = px.box(y, title=title, points=False)      # no outliers, whiskers go to min/max
    fig.update_layout(xaxis_title=elementtype, yaxis_title=None,
                      xaxis=dict(tickmode='array', tickangle=90, tickvals=x,
                                 ticktext=x.str.replace(elementtype + ': ', '')))

    y_max = y.max()
    y_q3 = y.quantile(0.75)
    y_med = y.median()
    y_q1 = y.quantile(0.25)
    y_min = y.min()
    diffmax = y_max.max() - y_min.min()
    # fig2 hides and replaces the hover information from fig
    fig2 = px.bar(x=x, y=y_max-y_min+0.1*diffmax, base=y_min-0.05*diffmax,
                  custom_data=[y_max, y_q3, y_med, y_q1, y_min])
    fig2.update_traces(opacity=0,
                       hovertemplate="<br>".join([elementtype + ": %{x}",
                            "max = %{customdata[0]}", "Q3 = %{customdata[1]}", "median = %{customdata[2]}",
                            "Q1 = %{customdata[3]}", "min = %{customdata[4]}"]))
    fig.add_traces(list(fig2.data))
    fig.update_layout(bargap=0)
    return fig


def compute_tree(df, df_len=None, df_diam=None):
    # https://python.plainenglish.io/create-a-network-graph-in-python-8829e0ec6741
    # https://stackoverflow.com/questions/57512155/how-to-draw-a-tree-more-beautifully-in-networkx

    lengths = diams = [None] * (len(df.index) + 1)      # used if df_len=None, df_diam=None

    neutrons = df.loc['n'].values
    neutrons = np.sort(np.unique(neutrons))
    neutrons = neutrons[neutrons > 0]
    neutrons = ['Bifurcation: '+str(neutron) for neutron in neutrons]

    G = nx.Graph()
    G.add_node('Central Plant', pos=(0, 1))         # Central Plant
    if len(neutrons) > 0:
        G.add_node('Bifurcation: 1', pos=(0, 0))    # first neutron

    # Branch 0
    protons = df.loc['k', '0'].values
    protons = protons[protons > 0]
    protons = ['Tee: ' + str(proton) for proton in protons]
    if df_len is not None:
        lengths = df_len['0']
    if df_diam is not None:
        diams = df_diam['0']
    if len(protons) > 0:
        for i in range(len(protons)):
            G.add_node(protons[i], pos=(0, 1 - (i + 1) / (len(protons) + 1)))
            if i == 0:
                G.add_edges_from([('Central Plant', protons[i])],
                                 label="Segment: 0,1", length=lengths[i], diam=diams[i])
            else:
                G.add_edges_from([(protons[i - 1], protons[i])],
                                 label="Segment: 0,{}".format(i+1), length=lengths[i], diam=diams[i])
        if 'Bifurcation: 1' in G.nodes:
            G.add_edges_from([(protons[-1], 'Bifurcation: 1')],
                             label="Segment: 0,{}".format(len(protons)+1),
                             length=lengths[len(protons)], diam=diams[len(protons)])
    elif 'Bifurcation: 1' in G.nodes:
        G.add_edges_from([('Central Plant', 'Bifurcation: 1')],
                         label="Segment: 0,1", length=lengths[0], diam=diams[0])

    # Branches below neutrons
    for n, neutron in enumerate(neutrons, 1):
        pos = nx.get_node_attributes(G, 'pos')
        current_pos = pos[neutron]
        factor = pow(0.5, (-current_pos[1]))
        for sign in [-1, 1]:    # left/right branch below neutron
            next_pos = (current_pos[0] + sign*factor, current_pos[1] - 1)
            protons = df.loc['k', str(sign*n)].values
            protons = protons[protons > 0]
            protons = ['Tee: ' + str(proton) for proton in protons]
            if df_len is not None:
                lengths = df_len[str(sign*n)]
            if df_diam is not None:
                diams = df_diam[str(sign*n)]
            for i in range(len(protons)):
                G.add_node(protons[i],
                           pos=(current_pos[0] + (i+1)/(len(protons)+1)*(next_pos[0] - current_pos[0]),
                                current_pos[1] + (i+1)/(len(protons)+1)*(next_pos[1] - current_pos[1])))
                if i == 0:
                    G.add_edges_from([(neutron, protons[i])],
                                     label="Segment: {},1".format(sign*n), length=lengths[i], diam=diams[i])
                else:
                    G.add_edges_from([(protons[i-1], protons[i])],
                                     label="Segment: {},{}".format(sign*n, i+1), length=lengths[i], diam=diams[i])
            neutron_below = df.loc['n', str(sign*n)].values[0]
            if neutron_below > 0:
                neutron_below = 'Bifurcation: ' + str(neutron_below)
                G.add_node(neutron_below, pos=next_pos)
                if len(protons) > 0:
                    G.add_edges_from([(protons[-1], neutron_below)],
                                     label="Segment: {},{}".format(sign*n, len(protons)+1),
                                     length=lengths[len(protons)], diam=diams[len(protons)])
                else:
                    G.add_edges_from([(neutron, neutron_below)],
                                     label="Segment: {},1".format(sign*n),
                                     length=lengths[len(protons)], diam=diams[len(protons)])

    # pos = nx.get_node_attributes(G, 'pos')
    # nx.draw(G, pos, with_labels=True)
    # plt.show()
    return G


def plot_tree(G):
    edge_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines',
        name='Segments')

    midpoint_trace = go.Scatter(
        x=[],
        y=[],
        text=[],
        mode='none',
        hoverinfo='text',
        # marker_size=0.5,
        # text=[0.45, 0.7, 0.34],
        # textposition='top center',
        # hovertemplate='weight: %{text}<extra></extra>',
        showlegend=False
    )

    for edge in G.edges(data=True):
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        edge_trace['x'] += tuple([x0, x1, None])
        edge_trace['y'] += tuple([y0, y1, None])
        midpoint_trace['x'] += tuple([(x0 + x1) / 2])
        midpoint_trace['y'] += tuple([(y0 + y1) / 2])
        edge_text = edge[2]['label']
        if edge[2]['length'] is not None:
            edge_text += '<br>Length = ' + str(edge[2]['length'])
        if edge[2]['diam'] is not None:
            edge_text += '<br>Diameter = ' + str(edge[2]['diam'])
        midpoint_trace['text'] += tuple([edge_text])

    fig = go.Figure()
    colors = fig.layout.template.layout.colorway  # get current color scheme

    central_trace = go.Scatter(
        x=[], y=[], text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            color=colors[2],
            size=20,
            symbol='circle'),
        name='Central Plant')

    neutron_trace = go.Scatter(
        x=[], y=[], text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            color=colors[0],
            size=15,
            symbol='circle'),
        name='Bifurcations')

    proton_trace = go.Scatter(
        x=[], y=[], text=[],
        mode='markers',
        hoverinfo='text',
        marker=dict(
            color=colors[1],
            size=10,
            symbol='circle'),
        name='Tees')

    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        if node == 'Central Plant':
            central_trace['x'] += tuple([x])
            central_trace['y'] += tuple([y])
            central_trace['text'] += tuple([node])
        elif node.startswith('Bifurcation'):
            neutron_trace['x'] += tuple([x])
            neutron_trace['y'] += tuple([y])
            neutron_trace['text'] += tuple([node])
        else:
            proton_trace['x'] += tuple([x])
            proton_trace['y'] += tuple([y])
            proton_trace['text'] += tuple([node])

    fig = go.Figure(data=[edge_trace, midpoint_trace, central_trace, neutron_trace, proton_trace],
                    layout=go.Layout(
                        title="Quantum Network",
                        titlefont=dict(size=16),
                        # showlegend=True,       should be called for each trace individually
                        legend=dict(
                            yanchor="top", y=1, xanchor="right", x=1, bgcolor='rgba(255,255,255,0.5)'),
                        hovermode='closest',
                        margin=dict(b=21, l=5, r=5, t=40),
                        # annotations=[dict(
                        #     text="Text Here",
                        #     showarrow=False,
                        #     xref="paper", yref="paper")],
                        xaxis=dict(visible=False, showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(visible=False, showgrid=False, zeroline=False, showticklabels=False)
                        )
                    )
    # fig.show()
    # plotly.offline.plot(fig)
    return fig


def get_treeinfo(G):
    neutrons = [node for node in G.nodes if node.startswith('Bifurcation')]
    protons = [node for node in G.nodes if node.startswith('Tee')]
    lengths = nx.get_edge_attributes(G, "length").values()
    if None not in lengths:
        total_length = sum(lengths)
    else:
        total_length = 'not specified'
    treeinfo = {"Number of Neutrons": len(neutrons),
                "Number of Protons": len(protons),
                "Number of Segments": len(G.edges),
                "Total Length": total_length}
    return treeinfo

'''--------------------------------------------------------------------------------------------------------'''

# interesting interactive solution using FigureWidget:
# https://github.com/plotly/plotly.py/issues/1445
def plot_widget(x, y, name=None, violin=False, colsize=0.9, spacing=0.01):
    # make y list of arrays
    if len(x) == len(y):
        y = [y]
    dim = len(y)
    if name is None:
        name = [None] * dim
    elif not isinstance(name, list):
        name = [name]
    if dim > len(name):
        name = name + [None] * (dim - len(name))
    elif dim < len(name):
        name = name[:len(y)]

    show_legend = [(s is not None) for s in name]

    fig = go.FigureWidget()
    colors = fig.layout.template.layout.colorway  # get current color scheme

    for i in range(dim):
        fig.add_scatter(x=x, y=y[i], mode='lines', name=name[i], line_color=colors[i], showlegend=show_legend[i])
        if violin:
            fig.add_violin(y=y[i], name=name[i], xaxis='x2', side='positive', points=False,
                           meanline=dict(visible=True), width=5, line_color=colors[i], showlegend=False)

    if violin:
        fig.layout = dict(xaxis=dict(domain=[0, colsize-spacing],
                                     showline=True, mirror=True),
                          yaxis=dict(showline=True, mirror='ticks'),
                          hovermode='closest',
                          xaxis2=dict(domain=[colsize+spacing, 1],
                                      showline=False, ticks='', showticklabels=False),
                          yaxis2=dict(showline=False, ticks='', showticklabels=False),
                          legend=dict(yanchor="top", y=1, xanchor="right", x=colsize-spacing,
                                      bgcolor='rgba(255,255,255,0.5)'),
                          )
    else:
        fig.layout = dict(xaxis=dict(domain=[0, 1], showline=True, mirror=True),
                          yaxis=dict(domain=[0, 1], showline=True, mirror='ticks'),
                          hovermode='closest',
                          legend=dict(yanchor="top", y=1, xanchor="right", x=1,
                                      bgcolor='rgba(255,255,255,0.5)'),
                          )

    def zoom(layout, xaxis_range, yaxis_range):
        for i in range(dim):
            idx = ((xaxis_range[0] <= x) & (x <= xaxis_range[1]) &
                   (yaxis_range[0] <= y[i]) & (y[i] <= yaxis_range[1]))
            # every second trace needs to be modified (violins)
            fig['data'][2*i+1]['y'] = y[i][idx]

    if violin:
        fig.layout.on_change(zoom, 'xaxis.range', 'yaxis.range')    # no syntax problem here!
        fig.update_traces(selector=dict(type="violin"),
                          hoverinfo='y', hoveron='kde')

    return fig


# create similar plot using subplots (use in case if WebGL3 not supported)
def plot_subplot(x, y, name=None, violin=False, colsize=0.9, spacing=0.01):
    # make y list of arrays
    if len(x) == len(y):
        y = [y]
    dim = len(y)
    if name is None:
        name = [None] * dim
    elif not isinstance(name, list):
        name = [name]
    if dim > len(name):
        name = name + [None] * (dim - len(name))
    elif dim < len(name):
        name = name[:len(y)]

    show_legend = [(s is not None) for s in name]

    if violin:
        colwidth = [colsize, 1 - colsize]
    else:
        colwidth = [1]

    fig = make_subplots(rows=1, cols=len(colwidth), shared_yaxes=True, column_widths=colwidth,
                        horizontal_spacing=2*spacing)
    colors = fig.layout.template.layout.colorway    # get current color scheme

    for i in range(dim):
        fig.add_trace(go.Scatter(x=x, y=y[i],
                                 mode='lines', name=name[i], line_color=colors[i], showlegend=show_legend[i]),
                      row=1, col=1)
        if violin:
            fig.add_trace(go.Violin(y=y[i], side='positive', line_color=colors[i], showlegend=False),
                          row=1, col=2)

    fig.update_xaxes(showline=True, mirror=True, row=1, col=1)
    fig.update_yaxes(showline=True, mirror='ticks', row=1, col=1)
    if violin:
        fig.update_traces(meanline_visible=True, width=5, row=1, col=2)
        fig.update_xaxes(showline=False, ticks='', showticklabels=False, row=1, col=2)
        fig.update_yaxes(showline=False, ticks='', showticklabels=False, row=1, col=2)

    fig.update_layout(legend=dict(yanchor="top", y=1, xanchor="right", x=colwidth[0]-0.02,
                                  bgcolor='rgba(255,255,255,0.5)'))
    return fig
