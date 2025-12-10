import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import sys
import io
sys.path.insert(0, 'C:/Users/cornelia.blanke/Documents/adv/Python')
import myplotting
from myplotting import *
from mypdf import *

import plotly.io as pio
import plotly.graph_objs as go

from shiny import ui, render, App, reactive
from shiny.types import FileInfo
from shinywidgets import output_widget, render_widget

from datetime import timedelta


import re
import numbers


properties = ["massflow", "speed", "mu", "mubar", "dP", "dP_m",
              "demand", "power", "dPower", "TIn", "TOut", "dT", "delta", "nubar", "epsilon",
              "material_cost", "engineering_cost", "operating_cost", "resources_cost"]

keys = ["Case Directory", "Case Name", "Meteo File", "Meteo Year", "Air Temp", "Ground Temp",
        "Fluid", "Boiler Type", "Boiler Power", "Boiler Price/kWh",
        "Supply Temperature", "Estimated Return Temperature",
        "Solve Mass Flow", "Solve Temperature", "Solver", "Update COP", "Solve Cost",
        "Save Results", "Save Frequency"]   # + properties (if needed)

app_ui = ui.page_fluid(
    # ui.h2("Playing with colormaps"),
    # ui.markdown("""
    #     This app is based on a [Matplotlib example][0] that displays 2D data
    #     with a user-adjustable colormap. We use a range slider to set the data
    #     range that is covered by the colormap.
    #
    #     [0]: https://matplotlib.org/3.5.3/gallery/userdemo/colormap_interactive_adjustment.html
    # """),
    ui.navset_tab(
        ui.nav("Overview",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    ui.input_file("file", "Choose Settings File", accept=[".txt"], multiple=False),
                    ui.row('Case Directory', style="padding-left:12px"),
                    ui.output_text_verbatim("casedir", placeholder=True),
                    ui.row('Case Name', style="padding-left:12px"),
                    ui.output_text_verbatim("casename", placeholder=True),
                    ui.row(
                        ui.column(6, ui.input_action_button("loadbtn", "Load Files")),
                        ui.column(6, ui.input_action_button("addbtn", "Add Files")),
                    ),
                    ui.panel_conditional("input.loadbtn || input.addbtn",
                                         ui.row('Select Case (Display/Download)',
                                                style="padding-left:12px; padding-top:20px"),
                                         ui.input_select("netcase", "", choices=[]),
                                         ui.download_button("download", "Download pdf"))
                ),
                ui.panel_main(
                    ui.panel_conditional("input.loadbtn || input.addbtn",
                        output_widget("treeplot"),
                        ui.row(
                            ui.column(6, ui.output_table("table")),
                            ui.column(6, ui.output_table("table2")),
                        )
                    )
                )
            )
        ),
        ui.nav("Plot 1",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    #ui.panel_conditional("input.loadbtn",
                        ui.input_checkbox_group("case", "Select Case", choices=[]),
                        ui.input_select("property", "Select Property", choices=[]),
                        ui.row(
                            ui.column(6, ui.input_select("elementtype", "Select Element Type", choices=[])),
                            ui.column(6, ui.output_ui("idControls")),
                        ),
                        ui.input_radio_buttons("plottype", "Select plot type",
                                               {"plot_over_time": "Time-dependent data",
                                                "plot_over_time_sorted": "Sorted data",
                                                "plot_over_time_both": "Time-dependent and sorted data",
                                                "plot_over_time_cumulated": "Cumulated data"}),
                        ui.input_checkbox("violin", "Show Kernel Density Estimate")
                    #)
                ),
                ui.panel_main(
                    ui.panel_conditional("input.loadbtn || input.addbtn",
                        output_widget("plot"),
                    )
                )
            )
        ),
        ui.nav("Plot 2",
            ui.layout_sidebar(
                ui.panel_sidebar(
                    #ui.panel_conditional("input.loadbtn",
                        ui.input_checkbox_group("case2", "Select Case", choices=[]),
                        ui.input_select("property2", "Select Property", choices=[]),
                        ui.input_select("elementtype2", "Select Element Type", choices=[]),
                        ui.input_radio_buttons("plottype2", "Select Plot Type",
                                               {"plot_elements": "Stem Plot",
                                                "plot_boxplot": "Box Plot"}),
                        ui.input_checkbox("sort", "Sort Values"),
                        ui.panel_conditional("input.sort",
                            ui.row(
                                ui.column(4, "Sort by"),
                                ui.column(4, ui.input_select("sortmethod", "",
                                                {"mean": 'Mean', "min": 'Min', "max": 'Max'})),
                                style="padding-left:24px"
                            ),
                            ui.row(
                                ui.column(4, "Show"),
                                ui.column(4, ui.input_select("top", "", {"1": 'Highest', "-1": 'Lowest'})),
                                ui.column(4, ui.input_numeric("topnumber", "", min=0, step=1, value=0)),
                                style="padding-left:24px"
                        ))
                        #)
               ),
               ui.panel_main(
                   ui.panel_conditional("input.loadbtn || input.addbtn",
                        output_widget("plot2"),
                    )
                )
            )
        )
    )
)


def server(input, output, session):
    @reactive.Calc
    def load_settings():
        if input.file() is None:
            return None
        f: list[FileInfo] = input.file()
        with open(f[0]["datapath"], 'r') as file:
            content = file.readlines()      # list of lines as strings
            dicts = dict(zip(['Author'] + keys, ['']*(len(keys)+1)))     # dict with empty strings
            tmpkeys = keys.copy()
        for i, line in enumerate(content):
            if line.startswith('#Author'):
                dicts["Author"] = content[i+1].strip()    # remove leading/trailing \s \t \n
            else:
                for key in tmpkeys:
                    if line.startswith(key):
                        dicts[key] = line.removeprefix(key).strip()
                        tmpkeys.remove(key)
                        if key == "Meteo Year":     # old version
                            dicts["Air Temp"] = dicts["Meteo Year"]
                            dicts["Ground Temp"] = dicts["Meteo Year"]
                            tmpkeys.remove("Air Temp")
                            tmpkeys.remove("Ground Temp")
                        break
        del dicts["Meteo Year"]
        return dicts

    @reactive.Calc
    def load_files():
        input.loadbtn()
        df_all = {}  # dictionary of dictionaries
        with reactive.isolate():
            content = load_settings()
            if content is not None:
                casename = content["Case Name"]
            else:
                casename = None
            if casename is not None:
                df = read_files()
                df_all[casename] = df
        return df_all

    @reactive.Calc
    def add_files():
        input.addbtn()
        df_all = load_files()
        with reactive.isolate():
            content = load_settings()
            if content is not None:
                casename = content["Case Name"]
            else:
                casename = None
            if casename is not None and casename not in df_all.keys():
                df = read_files()
                if len(df.keys()) > 0:
                    df_all[casename] = df
        return df_all

    def read_files():
        df = {"settings": load_settings()}  # dictionary for settings-dict and the dataframes
        if df["settings"] is not None:
            casedir = df["settings"]["Case Directory"]
            casename = df["settings"]["Case Name"]
            file = os.path.join(casedir, casename + ".csv")
            if os.path.exists(file):
                df["network"] = pd.read_csv(file, sep='\t', index_col=[0, 1])
                df["network"] = df["network"].rename(columns={'0.1': '0'})

            file = os.path.join(casedir, casename + "_SST" + ".csv")
            if os.path.exists(file):
                df["SST"] = pd.read_csv(file, sep='\t',
                                        usecols=['Substation nb', 'Substation Type', 'L_raccord (m)', 'DN_raccord'],
                                        index_col='Substation nb', encoding='latin-1')

            for prop in ["length", "diam"]:
                file = os.path.join(casedir, casename + "_" + prop + ".csv")
                if os.path.exists(file):
                    df[prop] = pd.read_csv(file, sep='\t', index_col=0)
            for prop in properties:
                file = os.path.join(casedir, casename + "_results_" + prop + ".csv")
                if os.path.exists(file):
                    df[prop] = pd.read_csv(file, sep='\t').set_index("Timestep")

            # speed estimation only for segments and Fluid=Water:
            DENSITY = 975.0     # hard-coded approximation !!
            CP = 4200.0         # hard-coded approximation !!
            if "massflow" in df.keys() and "length" in df.keys() and "diam" in df.keys() \
                    and "speed" not in df.keys() and df["settings"]["Fluid"] == "Water":
                segments = df["massflow"].columns[df["massflow"].columns.str.contains("Segment")]
                if len(segments) > 0:
                    indices = [re.split(r'[,\s]', segment)[1:] for segment in segments]
                    diams = np.array([df["diam"][idx[0]].iloc[int(idx[1])-1]/1000 for idx in indices])
                    df["speed"] = df["massflow"][segments] / (DENSITY * diams**2/4 * 3.1415)
            # dP_m for segments:
            if "dP" in df.keys() and "length" in df.keys() and "dP_m" not in df.keys():
                segments = df["dP"].columns[df["dP"].columns.str.contains("Segment")]
                if len(segments) > 0:
                    indices = [re.split(r'[,\s]', segment)[1:] for segment in segments]
                    lengths = [df["length"][idx[0]].iloc[int(idx[1]) - 1] for idx in indices]
                    df["dP_m"] = df["dP"][segments] / (2 * lengths)     # two pipes !
            # power estimation only for Fluid=Water:
            if "massflow" in df.keys() and "dT" in df.keys()\
                    and "power" not in df.keys() and df["settings"]["Fluid"] == "Water":
                elements = [element for element in df["massflow"].columns if element in df["dT"].columns]
                df["power"] = df["massflow"][elements] * CP * df["dT"][elements]
            # speed estimation for SST and only for Fluid=Water
            df_SST = pd.DataFrame()
            file = os.path.join(casedir, casename + "_SST" + ".csv")
            if "massflow" in df.keys() and df["settings"]["Fluid"] == "Water" and \
                ("speed" not in df.keys() or (df["speed"].columns.str.contains("Substation").any() == False)):
                if os.path.exists(file):
                    df_SST = pd.read_csv(file, sep='\t',
                                         usecols=['Substation nb', 'L_raccord (m)', 'DN_raccord'],
                                         index_col='Substation nb')
                    SST = df["massflow"].columns[df["massflow"].columns.str.contains("Substation")]
                    for sst in SST:
                        diam = df_SST.at[int(sst.split()[1]), 'DN_raccord'] / 1000
                        if diam > 0:
                            df["speed"][sst] = df["massflow"][sst] / (DENSITY * diam**2/4 * 3.1415)
            # dP_m for SST
            if "dP" in df.keys() and \
                ("dP_m" not in df.keys() or (df["dP_m"].columns.str.contains("Substation").any() == False)):
                if 'L_raccord (m)' not in df_SST.columns and os.path.exists(file):
                    df_SST = pd.read_csv(file, sep='\t',
                                         usecols=['Substation nb', 'L_raccord (m)'],
                                         index_col='Substation nb')
                if 'L_raccord (m)' in df_SST.columns:   # successfully read
                    SST = df["dP"].columns[df["dP"].columns.str.contains("Substation")]
                    for sst in SST:
                        length = df_SST.at[int(sst.split()[1]), 'L_raccord (m)']
                        if length > 0:
                            df["dP_m"][sst] = df["dP"][sst] / (2 * length)
        return df

    @reactive.Calc
    def tree():
        df = add_files()
        case = input.netcase()
        G = treeinfo = df_len = df_diam = None
        if case is not None and case in df.keys() and "network" in df[case].keys():
            if "length" in df[case].keys():
                df_len = df[case]["length"]
            if "diam" in df[case].keys():
                df_diam = df[case]["diam"]
            G = compute_tree(df[case]["network"], df_len=df_len, df_diam=df_diam)
            treeinfo = get_treeinfo(G)
        return G, treeinfo

    @output
    @render.ui
    @reactive.event(input.elementtype)
    def idControls():
        if input.elementtype() == 'Segment':
            return ui.TagList(ui.row(
                ui.column(6, ui.input_select("id1", "l", choices=[])),
                ui.column(6, ui.input_select("id2", "k", choices=[]))
            ))
        elif input.elementtype() == 'Branch':
            return ui.TagList(ui.row(
                ui.column(6, ui.input_select("id1", "l", choices=[]))
            ))
        elif input.elementtype() == 'Bifurcation':
            return ui.TagList(ui.row(
                ui.column(6, ui.input_select("id1", "n", choices=[]))
            ))
        elif (input.elementtype() == 'Total') | (input.elementtype() == 'Central Plant'):
            return None
        else:   # Tee, Substation
            return ui.TagList(ui.row(
                ui.column(6, ui.input_select("id1", "ID", choices=[]))
            ))

    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0 and input.property() in df[cases[0]].keys():
            elementtype = list(df[cases[0]][input.property()].columns.str.split(':').str[0].unique())
            if input.elementtype() in elementtype:
                selection = input.elementtype()
            else:
                selection = elementtype[0]
            ui.update_select(
                "elementtype",
                choices=elementtype,
                selected=selection,
            )

    @reactive.Effect()
    def _():
        df = add_files()
        cases = list(df.keys())
        elementtype = input.elementtype()
        if len(cases) > 0 and input.property() in df[cases[0]].keys() and elementtype is not None:
            columns = df[cases[0]][input.property()].columns
            elements = columns[columns.str.contains(elementtype)]
            if len(elements) > 0:
                indices = [re.split(r'[,\s]', element) for element in elements]
                if elementtype in ['Total', 'Central Plant']:
                # if len(indices) == 1:
                    ui.update_select(
                        "id1",
                        choices=[],
                    )
                    ui.update_select(
                        "id2",
                        choices=[],
                    )
                elif len(indices[0]) > 1:
                    identities = [idx[1] for idx in indices]
                    if input.id1() in identities:
                        selection = input.id1()
                    else:
                        selection = identities[0]
                    ui.update_select(
                        "id1",
                        choices=identities,
                        selected=selection,
                    )
                    if len(indices[0]) == 2:
                        ui.update_select(
                            "id2",
                            choices=[],
                        )
                    else:
                        identities = [idx[2] for idx in indices if idx[1] == input.id1()]
                        if len(identities) > 0:
                            if input.id2() in identities:
                                selection = input.id2()
                            else:
                                selection = identities[0]
                            ui.update_select(
                                "id2",
                                choices=identities,
                                selected=selection,
                            )
    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0 and input.property2() in df[cases[0]].keys():
            elementtype = list(df[cases[0]][input.property2()].columns.str.split(':').str[0].unique())
            if input.elementtype2() in elementtype:
                selection = input.elementtype2()
            else:
                selection = elementtype[0]
            ui.update_select(
                "elementtype2",
                choices=elementtype,
                selected=selection,
            )

    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0:
            props = list(df[cases[0]].keys())
            props = [prop for prop in props if prop in properties[:-4]]         # no costs
            if input.property() in props:
                selection = input.property()
            elif len(props) > 0:
                selection = props[0]
            else:
                selection = []
            ui.update_select(
                "property",
                choices=props,
                selected=selection,
            )

    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0:
            props = list(df[cases[0]].keys())
            props = [prop for prop in props if prop in properties]
            if input.property2() in props:
                selection = input.property2()
            elif len(props) > 0:
                selection = props[0]
            else:
                selection = []
            ui.update_select(
                "property2",
                choices=props,
                selected=selection,
            )

    @reactive.Effect()
    @reactive.event(input.sort)
    def _():
        if input.sort():
            cases = input.case2()
            elementtype = input.elementtype2()
            prop = input.property2()
            if len(cases) > 0 and elementtype is not None:
                df = add_files()
                if prop in df[cases[0]].keys():
                    allvalue = len(df[cases[0]][prop].columns[df[cases[0]][prop].columns.str.contains(elementtype)])
                    ui.update_numeric("topnumber", value=allvalue)
        else:
            ui.update_numeric("topnumber", value=0)

    @reactive.Effect()
    def _():
        topnumber = input.topnumber()
        if topnumber is not None and topnumber != round(topnumber):
            ui.update_numeric("topnumber", value=round(topnumber))

    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0:
            ui.update_select(
                "netcase",
                choices=cases,
                selected=cases[-1]
            )

    @reactive.Effect()
    def _():
        # df = load_files()
        df = add_files()
        cases = list(df.keys())
        if len(cases) > 0:
            ui.update_checkbox_group(
                "case",
                choices=cases,
                selected=cases[0]
            )
            ui.update_checkbox_group(
                "case2",
                choices=cases,
                selected=cases[0]
            )

    @output
    @render.text
    def casedir():
        content = load_settings()
        if content is None:
            return None
        else:
            return content["Case Directory"]

    @output
    @render.text
    def casename():
        content = load_settings()
        if content is None:
            return None
        else:
            return content["Case Name"]

    @output
    @render_widget
    def treeplot():
        fig = go.Figure(data=[], layout=go.Layout(plot_bgcolor='white',     # empty plot
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
        G, _ = tree()
        if G is not None:
            fig = plot_tree(G)
            fig.update_layout(title=str(input.netcase()) + " - Quantum Network")
        return fig

    @output
    @render.table
    def table():
        _, treeinfo = tree()
        case = input.netcase()
        df = add_files()
        nb_HX = nb_HP = ''
        if case in df.keys() and "SST" in df[case].keys():
            nb_HX = str(len(df[case]["SST"][df[case]["SST"]["Substation Type"] == "HX"]))
            nb_HP = str(len(df[case]["SST"][df[case]["SST"]["Substation Type"] == "HP"]))

        df_table = pd.DataFrame()
        if treeinfo is not None:
            d = {'Number of Substations<br> - Heat Exchanger<br> - Heat Pump':
                 str(treeinfo["Number of Protons"]) + "<br>" + nb_HX + "<br>" + nb_HP,
                 'Number of Branches': 2 * treeinfo["Number of Neutrons"] + 1,
                 'Number of Bifurcations': treeinfo["Number of Neutrons"],
                 'Number of Segments': treeinfo["Number of Segments"],
                 'Total Length': treeinfo["Total Length"]}
            df_table = pd.DataFrame(d.items())
            df_table = apply_styler(df_table, "Topology of Network")
        return df_table

    @output
    @render.table
    def table2():
        case = input.netcase()
        df_table = pd.DataFrame()
        df = add_files()
        if case is not None and case in df.keys():
            content = df[case]["settings"]
            d = {'Fluid': content["Fluid"],
                 'Boiler Type': content["Boiler Type"],
                 'Boiler Power': content["Boiler Power"],
                 'Supply Temperature': content["Supply Temperature"],
                 'Estimated Return Temperature': content["Estimated Return Temperature"]}
            df_table = pd.DataFrame(d.items())
            df_table = apply_styler(df_table, "Central Plant")
        return df_table

    def apply_styler(df, title):
        return df.style.set_table_attributes('class="dataframe shiny-table table w-auto"')\
            .hide(axis="index").hide(axis="columns") \
            .format(precision=2) \
            .set_table_styles([{
                    'selector': 'caption',
                    'props': 'caption-side: top; font-size:1.25em; font-weight: bold;'}])\
            .set_caption(title)

    @output
    @render_widget
    def plot():
        plt.close()     # close currently open figure
        df = load_files()
        fig = go.Figure(data=[], layout=go.Layout(plot_bgcolor='white',  # empty plot
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
        cases = input.case()
        element = input.elementtype()
        if element not in ["Total", "Central Plant"] and input.id1() is not None:
            element = element + ": " + input.id1()
            if element.startswith("Segment") and input.id2() is not None:
                element = element + "," + input.id2()

        # if len(cases) > 0 and input.property() in df[cases[0]].keys()\
        #         and element in df[cases[0]][input.property()].columns:
        #     plottype = input.plottype()
        #     func = getattr(myplotting, plottype)
        #     fig = func(input.property(), element, [df[input.case()[0]], df[input.case()[1]]], violin=input.violin())

        df_list = []
        for i in range(len(cases)):
            if input.property() in df[cases[i]].keys() \
                    and element in df[cases[i]][input.property()].columns:
                df_list.append(df[cases[i]])
        if len(df_list) > 0:
            plottype = input.plottype()
            func = getattr(myplotting, plottype)
            fig = func(input.property(), element, df_list, name=list(cases), violin=input.violin())
        return fig

    @output
    @render_widget
    def plot2():
        plt.close()  # close currently open figure
        df = load_files()
        fig = go.Figure(data=[], layout=go.Layout(plot_bgcolor='white',  # empty plot
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
        cases = input.case2()
        prop = input.property2()
        if len(cases) > 0 and prop in df[cases[0]].keys():
            elementtype = input.elementtype2()
            plottype = input.plottype2()
            nb = input.topnumber()
            if input.sort() and nb is not None:
                method = input.sortmethod()
                top = int(input.top()) * nb
            else:
                top = method = None
            func = getattr(myplotting, plottype)
            # avoid fig function being called too often
            if (nb != 0) == input.sort():
                fig = func(prop, elementtype, df[cases[0]], top=top, method=method)
        return fig

    @session.download(filename=lambda: f"Report {input.netcase()}.pdf")
    async def download():
        out = create_report()
        yield io.BytesIO(out).getvalue()

    def create_report():
        pio.templates.default = 'simple_white'
        df = add_files()
        case = input.netcase()
        content = df[case]["settings"]
        title = 'Report ' + case + '.pdf'
        author = content["Author"]

        pdf = PDF(title, author)

        '''Summarize inputs'''
        pdf.add_page()
        pdf.set_font('Helvetica')

        pdf.introtable()

        pdf.chapter_title('CASE DEFINITION')
        # pdf.chapter_text('This is a test.')
        # pdf.chapter_formula(r"$x^n + y^n = \frac{a}{b}$")

        G, treeinfo = tree()
        fig = plot_tree(G)
        fig.update_layout(title=None)
        pdf.chapter_image(fig, 'Quantum Network')

        nb_HX = nb_HP = ''
        if "SST" in df[case].keys():
            nb_HX = str(len(df[case]["SST"][df[case]["SST"]["Substation Type"] == "HX"]))
            nb_HP = str(len(df[case]["SST"][df[case]["SST"]["Substation Type"] == "HP"]))

        pdf.infotable('Size of Network',
                      (('Number of Substations\n - Heat Exchanger\n - Heat Pump',
                        str(treeinfo["Number of Protons"]) + "\n" + nb_HX + "\n" + nb_HP),
                       ('Number of Branches', 2 * treeinfo["Number of Neutrons"] + 1),
                       ('Number of Bifurcations', treeinfo["Number of Neutrons"]),
                       ('Number of Segments', treeinfo["Number of Segments"]),
                       ('Total Length [m]', treeinfo["Total Length"])
                       ))

        pdf.infotable('Central Plant',
                      (('Fluid', content["Fluid"]),
                       ('Boiler Type', content["Boiler Type"]),
                       ('Boiler Power [W]', content["Boiler Power"]),
                       ('Supply Temperature [°C]', content["Supply Temperature"]),
                       ('Estimated Return Temp. [°C]', content["Estimated Return Temperature"])))

        pdf.infotable('Weather Data',
                      (('Meteo File', content["Meteo File"]),
                       ('Air Temp', content["Air Temp"]),
                       ('Ground Temp', content["Ground Temp"])))

        df_meteo = read_meteo(content["Meteo File"], [content["Air Temp"], content["Ground Temp"]])
        pdf.chapter_image(plot_meteo(df_meteo), 'Weather Data')


        '''Solver Settings'''
        pdf.add_page()
        pdf.chapter_title('SIMULATION SETTINGS')
        pdf.infotable('Settings',
                      (('Solve Mass Flow', content["Solve Mass Flow"]),
                       ('Solve Temperature', content["Solve Temperature"]),
                       ('Temperature Solver', content["Solver"]),
                       ('Update COP', content["Update COP"]),
                       ('Solve Cost', content["Solve Cost"])))



        '''Results'''
        pdf.add_page()
        pdf.chapter_title('RESULTS')
        pdf.chapter_subtitle('Flow Distribution')

        def to_str(datalist, precision=2):
            strdata = []
            for datum in datalist:
                if isinstance(datum, numbers.Number):
                    strdata.append(str(round(datum, precision)))
                else:
                    strdata.append(datum)
            return strdata

        def identify_time(timestep):
            return str(timestep) + " (" + get_datetime(timestep - 1) + ")"

        data = ['not specified'] * 15
        if "massflow" in df[case].keys() and "Central Plant" in df[case]["massflow"].columns:
            df_CP = df[case]["massflow"]["Central Plant"]
            data[0] = df_CP.max()
            data[1] = identify_time(df_CP.idxmax())
            data[2] = df_CP.mean()
        if "speed" in df[case].keys():
            df_segment = df[case]["speed"].filter(like='Segment', axis=1)
            if df_segment.shape[1] > 0:
                max_per_segment = df_segment.max(axis=0)    # yields max of each segment
                mean_per_segment = df_segment.mean(axis=0)  # yields mean of each segment
                max_per_timestep = df_segment.max(axis=1)   # yields max of each timestep
                data[3] = max_per_segment.max()
                data[4] = max_per_segment.idxmax()
                data[5] = identify_time(max_per_timestep.idxmax())
                data[6] = mean_per_segment.max()
                data[7] = mean_per_segment.idxmax()
            df_SST = df[case]["speed"].filter(like='Substation', axis=1)
            if df_SST.shape[1] > 0:
                max_per_SST = df_SST.max(axis=0)            # yields max of each SST
                mean_per_SST = df_SST.mean(axis=0)          # yields mean of each SST
                max_per_timestep = df_SST.max(axis=1)       # yields max of each timestep
                data[8] = max_per_SST.max()
                data[9] = max_per_SST.idxmax()
                data[10] = identify_time(max_per_timestep.idxmax())
                data[11] = mean_per_SST.max()
                data[12] = mean_per_SST.idxmax()

        data = to_str(data)
        pdf.infotable('Flow Results',
                      (('Maximum Mass Flow from Central Plant [kg/s]\n - on Timestep (day/hour)',
                        data[0] + "\n" + data[1]),
                       ('Time-Averaged Mass Flow from Central Plant [kg/s]', data[2]),
                       ('Maximum Speed in Segments [m/s]\n - in Segment\n - on Timestep (day/hour)',
                        "\n".join(data[3:6])),
                       ('Maximum Time-Averaged Speed in Segments [m/s]\n - in Segment',
                        data[6] + "\n" + data[7]),
                       ('Maximum Speed in SST-Connections [m/s]\n - in Substation\n - on Timestep (day/hour)',
                        "\n".join(data[8:11])),
                       ('Maximum Time-Averaged Speed in SST-Connections [m/s]\n - in Substation',
                        data[11] + "\n" + data[12])
                       ),
                      col_widths=[0.7, 0.3])

        fig = plot_over_time_both("massflow", "Central Plant", df[case])
        fig.update_layout(title=None)
        pdf.chapter_image(fig, 'Mass Flow from Central Plant')

        pdf.chapter_text("""For a substation, the quantum number mubar is defined as the local mass flow
divided by the nominal mass flow per substation.""")
        pdf.chapter_formula(r"$\bar{\mu}(t) = \mu_\mathrm{SST}(t) = "
                            r"\frac{\mathrm{Mass\:Flow}\,(t)}"
                            r"{\frac{Nominal\:Mass\:Flow\:of\:CP}{Number\:of\:Substations}}$")
        pdf.chapter_text("""For a segment, mubar is then the arithmetic average of all substations alimented by that
segment.""")
        pdf.chapter_formula(r"$\bar{\mu}(t) ="
                            r"\frac{1}{Number\:of\:alimented\:SST}"
                            r"\sum_{\mathrm{SST}} \mu_\mathrm{SST} (t)$")

        fig = make_subplots(rows=1, cols=2, shared_yaxes=True, column_widths=[0.5, 0.5],
                            horizontal_spacing=0.1)
        fig1 = plot_boxplot("mubar", "Substation", df[case], top=10)
        fig2 = plot_boxplot("mubar", "Segment", df[case], top=10)
        fig.add_trace(fig1['data'][0], row=1, col=1)    # only boxplot (but not barplot)
        fig.add_trace(fig2['data'][0], row=1, col=2)
        fig1['layout']['xaxis'].pop('domain')           # take same layout but remove domain extend
        fig2['layout']['xaxis'].pop('domain')
        fig.update_xaxes(fig1['layout']['xaxis'], row=1, col=1)
        fig.update_xaxes(fig2['layout']['xaxis'], row=1, col=2)
        pdf.chapter_image(fig, 'mubar in Substations and Segments (Top 10 sorted by mean)')


        data = ['not specified'] * 15
        if "dP" in df[case].keys() and "Central Plant" in df[case]["dP"].columns:
            df_CP = df[case]["dP"]["Central Plant"]
            data[0] = df_CP.max()
            data[1] = identify_time(df_CP.idxmax())
            data[2] = df_CP.mean()
        if "dP_m" in df[case].keys():
            df_segment = df[case]["dP_m"].filter(like='Segment', axis=1)
            if df_segment.shape[1] > 0:
                max_per_segment = df_segment.max(axis=0)    # yields max of each segment
                mean_per_segment = df_segment.mean(axis=0)  # yields mean of each segment
                max_per_timestep = df_segment.max(axis=1)   # yields max of each timestep
                data[3] = max_per_segment.max()
                data[4] = max_per_segment.idxmax()
                data[5] = identify_time(max_per_timestep.idxmax())
                data[6] = mean_per_segment.max()
                data[7] = mean_per_segment.idxmax()
            df_SST = df[case]["dP_m"].filter(like='Substation', axis=1)
            if df_SST.shape[1] > 0:
                max_per_SST = df_SST.max(axis=0)            # yields max of each SST
                mean_per_SST = df_SST.mean(axis=0)          # yields mean of each SST
                max_per_timestep = df_SST.max(axis=1)       # yields max of each timestep
                data[8] = max_per_SST.max()
                data[9] = max_per_SST.idxmax()
                data[10] = identify_time(max_per_timestep.idxmax())
                data[11] = mean_per_SST.max()
                data[12] = mean_per_SST.idxmax()
        data = to_str(data)
        pdf.infotable('Pressure Loss Results',
                      (('Maximum Total Pressure Loss [Pa]\n - on Timestep (day/hour)',
                        data[0] + "\n" + data[1]),
                       ('Time-Averaged Total Pressure Loss [Pa]', data[2]),
                       ('Maximum Pressure Loss per meter [Pa/m]\n - in Segment\n - on Timestep (day/hour)',
                        "\n".join(data[3:6])),
                       ('Maximum Time-Averaged Pressure Loss per meter [Pa/m]\n - in Segment',
                        data[6] + "\n" + data[7]),
                       ('Maximum Pressure Loss per meter [Pa/m]\n - in Substation\n - on Timestep (day/hour)',
                        "\n".join(data[8:11])),
                       ('Maximum Time-Averaged Pressure Loss per meter [Pa/m]\n - in Substation',
                        data[11] + "\n" + data[12]),
                       ),
                      col_widths=[0.7, 0.3])


        pdf.chapter_subtitle('Energy Consumption')
        df_CP = pd.Series()
        df_SST = pd.DataFrame()
        data = ['not specified'] * 10
        if "power" in df[case].keys():
            if "Central Plant" in df[case]["power"].columns:
                df_CP = df[case]["power"]["Central Plant"]
                data[0] = df_CP.sum() /1e6
                data[5] = df_CP.max() / 1000
                data[6] = identify_time(df_CP.idxmax())
                data[7] = df_CP.mean() / 1000
            df_SST = df[case]["power"].filter(like='Substation', axis=1).copy()
            if df_SST.shape[1] > 0:
                df_SST["sum"] = df_SST.sum(axis=1)  # sum of all SST
                data[1] = df_SST["sum"].sum() / 1e6
            if isinstance(data[0], numbers.Number) and isinstance(data[1], numbers.Number):
                data[2] = (data[0] - data[1])
        if "demand" in df[case].keys():
            df_demand = df[case]["demand"].filter(like='Substation', axis=1)
            if df_demand.shape[1] > 0:
                data[3] = df_demand.sum().sum() / 1e6
            if isinstance(data[1], numbers.Number) and isinstance(data[3], numbers.Number):
                data[4] = (data[1] - data[3])

        data = to_str(data)
        pdf.infotable('Energy Results',
                      (('Total Energy provided by Central Plant [MWh]', data[0]),
                       ('Total Energy in Substations [MWh]', data[1]),
                       ('Total Energy Loss in Network [MWh]', data[2]),
                       ('Total Energy Demand from Consumers [MWh]', data[3]),
                       ('Total Energy Loss in Substations [MWh]', data[4]),
                       ('Maximum Provided Power [kW]\n - on Timestep (day/hour)\n - Time-Averaged Provided Power [kW]',
                        "\n".join(data[5:8]))),
                      col_widths=[0.7, 0.3])

        fig = plot_over_time_both("power", "Central Plant", df[case])
        fig.update_layout(title=None)
        pdf.chapter_image(fig, 'Power from Central Plant')

        fig = make_subplots(specs=[[{"secondary_y": True}]])
        x = df_CP.index
        if len(x) > 0:
            y_CP = df_CP.sort_values(ascending=False)
            fig.add_scatter(x=x, y=y_CP.values,
                            mode='lines', name='Central Plant', secondary_y=False, showlegend=True)
            if "sum" in df_SST.columns:
                y_SST = df_SST["sum"].loc[y_CP.index]   # use ordering of y_CP
                y_diff = y_CP.values - y_SST.values
                y_perc = (y_diff / np.absolute(y_CP.values)) * 100
                y_perc = np.where(np.absolute(y_perc) > 100, np.nan, y_perc)   # clip to -100...100
                fig.add_scatter(x=x, y=y_SST.values,
                                mode='lines', name='Substations', secondary_y=False, showlegend=True)
                fig.add_scatter(x=x, y=y_diff,
                                mode='lines', name='Losses CP-SST', secondary_y=False, showlegend=True)
                fig.add_scatter(x=x, y=y_perc,
                                mode='lines', name='Losses Percentage', secondary_y=True, showlegend=True)
        fig.update_xaxes(title_text="Time [h]", zeroline=True)
        fig.update_yaxes(title_text="Power", rangemode="tozero", secondary_y=False)
        fig.update_yaxes(title_text="Percentage [%]", rangemode="tozero", secondary_y=True)
        fig.update_layout(title=None, xaxis=dict(showline=True, mirror=True))
        pdf.chapter_image(fig, 'Power Comparison')



        pdf.chapter_subtitle('Costs')
        data = ['not specified'] * 16
        count = 0
        for prop in ["material_cost", "engineering_cost", "operating_cost", "resources_cost"]:
            if prop in df[case].keys():
                if "Total" in df[case][prop].columns:
                    data[0 + count] = df[case][prop]["Total"].iloc[0] / 1000
                if "Central Plant" in df[case][prop].columns:
                    data[1 + count] = df[case][prop]["Central Plant"].iloc[0] / 1000
                df_SST = df[case][prop].filter(like='Substation', axis=1)
                df_segment = df[case][prop].filter(like='Segment', axis=1)
                if df_SST.shape[1] > 0:
                    data[2 + count] = df_SST.sum(axis=1).iloc[0] / 1000
                if df_segment.shape[1] > 0:
                    data[3 + count] = df_segment.sum(axis=1).iloc[0] / 1000
            count += 4


        pdf.infotable('Costs Results',
                      (('in kCHF', 'Total', 'Central Plant', 'Substations', 'Segments'),
                       ('Material Cost', data[0], data[1], data[2], data[3]),
                       ('Engineering Cost', data[4], data[5], data[6], data[7]),
                       ('Operation per Year', data[8], data[9], data[10], data[11]),
                       ('Resources per Year', data[12], data[13], data[14], data[15])),
                      col_widths=[0.24, 0.19, 0.19, 0.19, 0.19])



        # reset default
        pio.templates.default = 'plotly'
        return pdf.output()

    def get_datetime(hour):
        #return format(pd.Timestamp(year) + timedelta(hours=float(hour)), '%d.%m.%Y %H:%M')
        return format(pd.Timestamp(1901) + timedelta(hours=float(hour)), '%d.%m. %H:%M')     #TODO leap year

app = App(app_ui, server)
