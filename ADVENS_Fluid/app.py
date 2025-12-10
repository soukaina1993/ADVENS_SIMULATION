import subprocess, math
import sys, os
import plotly.graph_objs as go
from shiny import App, Inputs, Outputs, Session, render, ui, reactive
from shinywidgets import output_widget, render_widget

if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    base_path = sys._MEIPASS
    os.environ['MY_BASE_PATH'] = base_path
else:
    base_path = ""


prop_keys = ["P", "T", "v", "h", "s", "u", "x"]
crit_keys = ["Pc", "Tc", "vc", "hc", "sc", "uc"]
var = dict(zip(prop_keys+crit_keys, range(len(prop_keys+crit_keys))))
prop_values = ["Pressure", "Temperature", "Specific Volume", "Enthalpy",
               "Entropy", "Inner Energy", "Vapour Quality"]
prop_units = ["[bar]", "[°C]", "[m3/kg]", "[kJ/kg]", "[kJ/kg.K]", "[kJ/kg]", "[-]"]
#prop_units = ["[Pa]", "[K]", "[m3/kg]", "[J/kg]", "[J/kg.K]", "[J/kg]", "[-]"]
prop_values_units = [m+' '+n for m, n in zip(prop_values, prop_units)]
prop_keys_units = [m+' '+n for m, n in zip(prop_keys, prop_units)]
props = dict(zip(prop_keys, prop_values_units))


ADVENS_Fluid = os.path.join(base_path, 'ADVENS_Fluid.exe')

PRFluidList = subprocess.run([ADVENS_Fluid, 'FluidsList', 'PRFluid'],
                             capture_output=True, text=True, close_fds=True).stdout.split(',')
PRFluidList = [fluid for fluid in PRFluidList if 'H2O' not in fluid]

CPFluidList = subprocess.run([ADVENS_Fluid, 'FluidsList', 'CoolProp'],
                             capture_output=True, text=True, close_fds=True).stdout.split(',')


# conversion to positive float, only if no-NaN
def to_float(instring):
    if instring == 'nan':
        raise ValueError
    else:
        value = float(instring)     # may rise ValueError
        if value < 0 and value != -1:
            raise ValueError
            #return value
        else:
            return value


def to_SI_units(prop, value):
    if prop == "P":
        return value * 1.0e5
    elif prop == "T":
        return value + 273.15
    elif prop == "h" or prop == "s" or prop == "u":
        return value * 1000.0
    else:
        return value

def from_SI_units(prop, value):
    if prop == "P":
        return value / 1.0e5
    elif prop == "T":
        return value - 273.15
    elif prop == "h" or prop == "s" or prop == "u":
        return value / 1000.0
    else:
        return value

# scale a list of values (for plotting)
def pscale(p_list):
    return [p/1.0e5 for p in p_list]
def hscale(h_list):         # for h and s
    return [h/1.0e3 for h in h_list]
def tscale(t_list):
    return [t-273.15 for t in t_list]



app_ui = ui.page_sidebar(
    ui.sidebar(
        ui.input_select("method", "Fluid Library",
                        {"PRFluid": "Peng Robinson", "MWater": "MWater", "CoolProp": "CoolProp"},
                        selected="CoolProp"),
        ui.input_select("fluid", "Fluid", []),
        ui.input_select("inprop1", "Property 1", {"P": "Pressure [bar]", "T": "Temperature [°C]"}),
        ui.input_numeric("invalue1", "Value:", 0),
        ui.input_select("inprop2", "Property 2", []),
        ui.input_numeric("invalue2", "Value:", 0),
        ui.input_action_button("compute", "Compute!"),

    ),
    ui.card(
        ui.navset_card_tab(*[
            ui.nav_panel("P-v", output_widget("Pv")),
            ui.nav_panel("P-h", output_widget("Ph")),
            ui.nav_panel("T-s", output_widget("Ts")),
            ui.nav_panel("h-s", output_widget("hs")),
        ]),
        ui.output_ui("html_table"),
    )
)


def server(input: Inputs, output: Outputs, session: Session):
    @reactive.calc
    @reactive.event(input.compute, ignore_none=True)
    def results():
        invalue1 = str(to_SI_units(input.inprop1(), input.invalue1()))
        invalue2 = str(to_SI_units(input.inprop2(), input.invalue2()))
        result = subprocess.run([ADVENS_Fluid,
                                 input.inprop1(), invalue1,
                                 input.inprop2(), invalue2,
                                 input.fluid(), input.method()],
                                capture_output=True, text=True, close_fds=True)
        if len(result.stderr) > 0:
            resultlist = ['nan']*7
            for i in range(len(prop_keys)):
                if prop_keys[i] == input.inprop1():
                    resultlist[i] = invalue1
                elif prop_keys[i] == input.inprop2():
                    resultlist[i] = invalue2
            resultlist += result.stdout.split(' ')[-6:]
            return resultlist
        else:
            return result.stdout[:-1].split(' ')        # remove '\n'

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def saturation_curve():
        resultlist = results()
        try:
            P = [to_float(resultlist[var['P']])]
        except ValueError:
            P = [float(resultlist[var['Pc']])]
        Pc = [float(resultlist[var['Pc']])]
        Tc = [float(resultlist[var['Tc']])]
        vc = [float(resultlist[var['vc']])]
        hc = [float(resultlist[var['hc']])]
        sc = [float(resultlist[var['sc']])]

        Pplot = [Pc[0]]
        Tplot = [Tc[0]]
        vliq = [vc[0]]
        vgas = [vc[0]]
        hliq = [hc[0]]
        hgas = [hc[0]]
        sliq = [sc[0]]
        sgas = [sc[0]]

        steps = 20
        logPmin = math.log(min(Pc[0] / 1000, P[0] / 10))
        logPc = math.log(Pc[0])

        # create list of pressures for evaluation
        Pevals = ','.join([str(math.exp(logPmin + (logPc - logPmin) * (steps - 1 - i) / steps)) for i in range(steps)])

        # compute all pressures for x=0 and x=1
        res = subprocess.run([ADVENS_Fluid,
                              'P', Pevals, 'x', '0,1', input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        for i in range(len(res)):
            try:
                if i % 2 == 0:                  # x = 0
                    Pplot.append(to_float(res[i][var['P']]))
                    Tplot.append(to_float(res[i][var['T']]))
                    vliq.append(to_float(res[i][var['v']]))
                    hliq.append(to_float(res[i][var['h']]))
                    sliq.append(to_float(res[i][var['s']]))
                else:                           # x = 1
                    vgas.append(to_float(res[i][var['v']]))
                    hgas.append(to_float(res[i][var['h']]))
                    sgas.append(to_float(res[i][var['s']]))
            except ValueError:
                break
        return Pplot, Tplot, vliq, vgas, hliq, hgas, sliq, sgas

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isobar_curve():
        resultlist = results()
        sreturn = []
        T = []
        h = []
        try:
            strP = resultlist[var['P']]
            to_float(strP)      # check if P is a float
            s = to_float(resultlist[var['s']])
            Tc = to_float(resultlist[var['Tc']])
        except ValueError:
            return sreturn, T, h

        # TODO P > Pc
        res = subprocess.run([ADVENS_Fluid,
                              'P', strP, 'x', '0,1', input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        sl = to_float(res[0][var['s']])
        Tl = to_float(res[0][var['T']])
        hl = to_float(res[0][var['h']])

        sg = to_float(res[1][var['s']])
        hg = to_float(res[1][var['h']])

        _, _, _, _, _, _, sliq, sgas = saturation_curve()
        smin = min(min(sliq), s)
        smax = max(max(sgas), s)

        steps = 10
        sevallist = [str(smin + (sl - smin) * i/steps) for i in range(steps)]               # liquid
        sevallist.extend([str(sg + (smax - sg) * i / steps) for i in range(1, steps + 1)])  # gas
        sevals = ','.join(sevallist)

        res = subprocess.run([ADVENS_Fluid,
                              'P', strP, 's', sevals, input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        for i in range(2*steps):
            if i == steps:
                sreturn.append(sl)
                T.append(Tl)
                h.append(hl)
                sreturn.append(sg)
                T.append(Tl)
                h.append(hg)
            try:
                Tact = to_float(res[i][var['T']])
                if Tact < Tc:
                    sreturn.append(to_float(res[i][var['s']]))
                    T.append(Tact)
                    h.append(to_float(res[i][var['h']]))
                else:
                    continue
            except ValueError:
                continue

        return sreturn, T, h

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isotherm_curve():
        resultlist = results()
        Preturn = []
        v = []
        s = []
        h = []
        try:
            strT = resultlist[var['T']]
            to_float(strT)
            Pc = to_float(resultlist[var['Pc']])
        except ValueError:
            return Preturn, v, s, h

        # TODO P > Pc
        res = subprocess.run([ADVENS_Fluid,
                              'T', strT, 'x', '0,1', input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists
        sl = to_float(res[0][var['s']])
        Pl = to_float(res[0][var['P']])
        hl = to_float(res[0][var['h']])
        vl = to_float(res[0][var['v']])

        sg = to_float(res[1][var['s']])
        hg = to_float(res[1][var['h']])
        vg = to_float(res[1][var['v']])

        logPmax = math.log(Pc * 0.99)
        logPmin = math.log(min(Pc / 1000, Pl / 10))
        logPl = math.log(Pl)

        steps = 10
        Pevallist = [str(math.exp(logPmax - i/steps * (logPmax - logPl))) for i in range(steps)]
        Pevallist.extend([str(math.exp(logPl - (logPl - logPmin) * i / steps)) for i in range(1, steps+1)])
        Pevals = ','.join(Pevallist)
        res = subprocess.run([ADVENS_Fluid,
                              'T', strT, 'P', Pevals, input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        for i in range(2*steps):
            if i == steps:
                Preturn.append(Pl)
                v.append(vl)
                s.append(sl)
                h.append(hl)
                Preturn.append(Pl)
                v.append(vg)
                s.append(sg)
                h.append(hg)
            try:
                Preturn.append(to_float(res[i][var['P']]))
                v.append(to_float(res[i][var['v']]))
                s.append(to_float(res[i][var['s']]))
                h.append(to_float(res[i][var['h']]))
            except ValueError:
                continue

        return Preturn, v, s, h

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isochore_curve():
        resultlist = results()
        Preturn = []
        Treturn = []
        s = []
        h = []
        try:
            strv = resultlist[var['v']]
            v = to_float(strv)
            P = to_float(resultlist[var['P']])
            Pc = to_float(resultlist[var['Pc']])
            Tc = to_float(resultlist[var['Tc']])
        except ValueError:
            return Preturn, Treturn, s, h

        # TODO P > Pc
        Tevallist = []
        _, Tsat, _, vgas, _, _, _, _ = saturation_curve()

        steps = 10
        Tmin = min(Tsat)
        Tmax = 0.99*Tc

        Tgas = 1e16
        for i in range(1, len(vgas)):
            if vgas[i] > v:
                factor = (v - vgas[i-1]) / (vgas[i] - vgas[i-1])
                Tgas = Tsat[i-1] - factor * (Tsat[i-1] - Tsat[i])
                break

        notfound = True
        for i in range(steps+1):
            Tact = Tmin + (Tmax - Tmin) * i / steps
            if notfound and Tact > Tgas:
                Tevallist.append(str(Tgas))
                notfound = False
            Tevallist.append(str(Tact))
        Tevals = ','.join(Tevallist)

        res = subprocess.run([ADVENS_Fluid,
                              'v', strv, 'T', Tevals, input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        notfound = True
        for result in res:
            Pact = to_float(result[var['P']])
            if notfound and Pact > P:    # append current results to list
                notfound = False
                Treturn.append(to_float(resultlist[var['T']]))
                Preturn.append(P)
                s.append(to_float(resultlist[var['s']]))
                h.append(to_float(resultlist[var['h']]))
            if Pact < Pc:
                try:
                    Treturn.append(to_float(result[var['T']]))
                    Preturn.append(Pact)
                    s.append(to_float(result[var['s']]))
                    h.append(to_float(result[var['h']]))
                except ValueError:
                    break
            elif len(Preturn) > 0:   # Pact >= Pc --> reduce Pact and compute final value
                Pact = (Preturn[-1] + Pc) / 2
                res = subprocess.run([ADVENS_Fluid,
                                      'v', strv, 'P', str(Pact), input.fluid(), input.method()],
                                     capture_output=True, text=True, close_fds=True).stdout[:-1].split(' ')
                try:
                    Treturn.append(to_float(res[var['T']]))
                    Preturn.append(Pact)
                    s.append(to_float(res[var['s']]))
                    h.append(to_float(res[var['h']]))
                except ValueError:
                    break
                break
        return Preturn, Treturn, s, h

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isentrope_curve():
        resultlist = results()
        Preturn = []
        T = []
        v = []
        h = []
        try:
            strs = resultlist[var['s']]
            to_float(strs)
            P = to_float(resultlist[var['P']])
            Tc = to_float(resultlist[var['Tc']])
        except ValueError:
            return Preturn, T, v, h

        # TODO P > Pc
        Psat, _, _, _, _, _, _, _ = saturation_curve()
        Psat[0] *= 0.99  # replace Pc with 0.99*Pc
        Pevals = ','.join([str(p) for p in Psat])

        res = subprocess.run([ADVENS_Fluid,
                              's', strs, 'P', Pevals, input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        notfound = True
        for result in res:
            try:
                Pact = to_float(result[var['P']])
                Tact = to_float(result[var['T']])
                if Tact < Tc:
                    try:
                        if notfound and P > Pact:
                            T.append(to_float(resultlist[var['T']]))
                            Preturn.append(P)
                            v.append(to_float(resultlist[var['v']]))
                            h.append(to_float(resultlist[var['h']]))
                            notfound = False
                        Preturn.append(to_float(result[var['P']]))
                        T.append(Tact)
                        v.append(to_float(result[var['v']]))
                        h.append(to_float(result[var['h']]))
                    except ValueError:
                        continue
            except ValueError:
                continue
        return Preturn, T, v, h

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isenthalpe_curve():
        resultlist = results()
        Preturn = []
        T = []
        v = []
        s = []
        try:
            strh = resultlist[var['h']]
            to_float(strh)
            P = to_float(resultlist[var['P']])
            Tc = to_float(resultlist[var['Tc']])
        except ValueError:
            return Preturn, T, v, s

        # TODO P > Pc
        Psat, _, _, _, _, _, _, _ = saturation_curve()
        Psat[0] *= 0.99  # replace Pc with 0.99*Pc
        Pevals = ','.join([str(p) for p in Psat])

        res = subprocess.run([ADVENS_Fluid,
                              'h', strh, 'P', Pevals, input.fluid(), input.method()],
                             capture_output=True, text=True, close_fds=True).stdout.split('\n')
        res = [s.split(' ') for s in res[:-1]]  # list of result lists

        notfound = True
        for result in res:
            try:
                Pact = to_float(result[var['P']])
                Tact = to_float(result[var['T']])
                if Tact < Tc:
                    try:
                        if notfound and P > Pact:
                            T.append(to_float(resultlist[var['T']]))
                            Preturn.append(P)
                            s.append(to_float(resultlist[var['s']]))
                            v.append(to_float(resultlist[var['v']]))
                            notfound = False
                        Preturn.append(Pact)
                        T.append(Tact)
                        s.append(to_float(result[var['s']]))
                        v.append(to_float(result[var['v']]))
                    except ValueError:
                        continue
                else:
                    continue
            except ValueError:
                continue
        return Preturn, T, v, s

    @reactive.calc
    @reactive.event(results, ignore_none=True)
    def isovapour_curve():
        resultlist = results()
        P = []
        T = []
        v = []
        s = []
        h = []
        try:
            x = to_float(resultlist[var['x']])
        except ValueError:
            return P, T, v, s, h
        if x < 0 or x > 1:
            return P, T, v, s, h

        P, T, vliq, vgas, hliq, hgas, sliq, sgas = saturation_curve()

        v = [vliq[i] + x * (vgas[i] - vliq[i]) for i in range(len(P))]
        s = [sliq[i] + x * (sgas[i] - sliq[i]) for i in range(len(P))]
        h = [hliq[i] + x * (hgas[i] - hliq[i]) for i in range(len(P))]
        return P, T, v, s, h


    @render.text
    @reactive.event(results, ignore_none=True)
    def html_table() -> str:
        resultlist = results().copy()       # do not override results!!
        if 'nan' in resultlist:
            html = '<b style="color:red">Error in Computation</b>'
        else:
            html = ''

        # conversion to non-SI units
        for i in range(len(resultlist)):
            resultlist[i] = str(from_SI_units(prop_keys[i % len(prop_keys)], to_float(resultlist[i])))
        resultlist += ['']

        html += '<table><tr><th style="padding-left:10px; padding-right:30px; border:1px solid black">Property</th>' + \
                '<th style="padding-left:10px; padding-right:30px; border:1px solid black">Unit</th>' + \
                '<th style="padding-left:10px; padding-right:30px; border:1px solid black">Value</th>' + \
                '<th style="padding-left:10px; padding-right:30px; border:1px solid black">Critical Point</th></tr>'
        for i in range(len(prop_values)):
            if prop_keys[i] == input.inprop1() or prop_keys[i] == input.inprop2():
                bgcolor = "lightgrey"
            else:
                bgcolor = "white"
            html += '<tr> <td style="padding-left:10px; padding-right:30px; border:1px solid black">' + \
                    prop_values[i] + '</td>' + \
                    '<td style="padding-left:10px; padding-right:30px; border:1px solid black" bgcolor=>' +\
                    prop_units[i] + '</td>' + \
                    '<td style="padding-left:10px; padding-right:30px; border:1px solid black" bgcolor=' + \
                    bgcolor + '>' + resultlist[i] + '</td>' + \
                    '<td style="padding-left:10px; padding-right:30px; border:1px solid black">' + \
                    resultlist[i+len(prop_values)] + '</td> </tr>'
        html += '</table>'
        return html


    @output
    @render_widget
    @reactive.event(results, ignore_none=True)
    def Pv() -> object:
        resultlist = results()
        try:
            v = [to_float(resultlist[var['v']])]
            P = [to_float(resultlist[var['P']])]
        except ValueError:
            v = []
            P = []
        vc = [float(resultlist[var['vc']])]
        Pc = [float(resultlist[var['Pc']])]

        Pplot, _, vliq, vgas,  _, _, _, _ = saturation_curve()
        Pisotherm, visotherm, _, _ = isotherm_curve()
        Pisentr, _, visentr, _ = isentrope_curve()
        Pisenth, _, visenth, _ = isenthalpe_curve()
        Pisovap, _, visovap, _, _ = isovapour_curve()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=vliq, y=pscale(Pplot), mode='lines', name='sat-liquid', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=vgas, y=pscale(Pplot), mode='lines', name='sat-gas', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=vc, y=pscale(Pc), mode='markers', name='critical point',
                                 marker=dict(size=10, color='black')))
        fig.add_trace(go.Scatter(x=v, y=pscale(P), mode='markers', name='result',
                                 marker=dict(size=15, color='magenta', symbol='star')))
        fig.add_trace(go.Scatter(x=[0, max(vgas)], y=pscale([P[0], P[0]]),
                                 mode='lines', name='isobar', line=dict(color='green'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=visotherm, y=pscale(Pisotherm), mode='lines', name='isotherm',
                                 line=dict(color='orange'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=[v[0], v[0]], y=pscale([min(Pplot), max(Pplot)]),
                                 mode='lines', name='isochore', line=dict(color='gold'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=visentr, y=pscale(Pisentr), mode='lines', name='isentrope',
                                 line=dict(color='aqua'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=visenth, y=pscale(Pisenth), mode='lines', name='isenthalpe',
                                 line=dict(color='lime'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=visovap, y=pscale(Pisovap), mode='lines', name='isovapour',
                                 line=dict(color='grey'), visible="legendonly"))
        fig.update_xaxes(title_text="v [m3/kg]", type="log")
        fig.update_yaxes(title_text="P [bar]", type="log")
        return fig

    @output
    @render_widget
    @reactive.event(results, ignore_none=True)
    def Ph() -> object:
        resultlist = results()
        try:
            h = [to_float(resultlist[var['h']])]
            P = [to_float(resultlist[var['P']])]
        except ValueError:
            h = []
            P = []
        hc = [float(resultlist[var['hc']])]
        Pc = [float(resultlist[var['Pc']])]

        Pplot, _, _, _, hliq, hgas,  _, _ = saturation_curve()
        Pisotherm, _, _, hisotherm = isotherm_curve()
        Pisochor, _, _, hisochor = isochore_curve()
        Pisentr, _, _, hisentr = isentrope_curve()
        Pisovap, _, _, _, hisovap = isovapour_curve()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=hscale(hliq), y=pscale(Pplot),
                                 mode='lines', name='sat-liquid', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=hscale(hgas), y=pscale(Pplot),
                                 mode='lines', name='sat-gas', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=hscale(hc), y=pscale(Pc), mode='markers', name='critical point',
                                 marker=dict(size=10, color='black')))
        fig.add_trace(go.Scatter(x=hscale(h), y=pscale(P), mode='markers', name='result',
                                 marker=dict(size=15, color='magenta', symbol='star')))
        fig.add_trace(go.Scatter(x=hscale([min(min(hliq), h[0]), max(max(hgas), h[0])]), y=pscale([P[0], P[0]]),
                                 mode='lines', name='isobar', line=dict(color='green'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(hisotherm), y=pscale(Pisotherm), mode='lines', name='isotherm',
                                 line=dict(color='orange'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(hisochor), y=pscale(Pisochor), mode='lines', name='isochore',
                                 line=dict(color='gold'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(hisentr), y=pscale(Pisentr), mode='lines', name='isentrope',
                                 line=dict(color='aqua'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale([h[0], h[0]]), y=pscale([min(Pplot), max(Pplot)]),
                                 mode='lines', name='isenthalpe', line=dict(color='lime'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(hisovap), y=pscale(Pisovap), mode='lines', name='isovapour',
                                 line=dict(color='grey'), visible="legendonly"))
        fig.update_xaxes(title_text="h [kJ/kg]")
        fig.update_yaxes(title_text="P [bar]", type="log")
        return fig

    @output
    @render_widget
    @reactive.event(results, ignore_none=True)
    def Ts() -> object:
        resultlist = results()
        try:
            s = [to_float(resultlist[var['s']])]
            T = [to_float(resultlist[var['T']])]
        except ValueError:
            s = []
            T = []
        sc = [float(resultlist[var['sc']])]
        Tc = [float(resultlist[var['Tc']])]

        _, Tplot, _, _, _, _, sliq, sgas = saturation_curve()
        sisobar, Tisobar, _ = isobar_curve()
        _, Tisochor, sisochor, _ = isochore_curve()
        _, Tisenth, _, sisenth = isenthalpe_curve()
        _, Tisovap, _, sisovap, _ = isovapour_curve()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=hscale(sliq), y=tscale(Tplot),
                                 mode='lines', name='sat-liquid', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=hscale(sgas), y=tscale(Tplot),
                                 mode='lines', name='sat-gas', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=hscale(sc), y=tscale(Tc), mode='markers', name='critical point',
                                 marker=dict(size=10, color='black')))
        fig.add_trace(go.Scatter(x=hscale(s), y=tscale(T), mode='markers', name='result',
                                 marker=dict(size=15, color='magenta', symbol='star')))
        fig.add_trace(go.Scatter(x=hscale(sisobar), y=tscale(Tisobar), mode='lines', name='isobar',
                                 line=dict(color='green'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale([min(min(sliq), s[0]), max(max(sgas), s[0])]), y=tscale([T[0], T[0]]),
                                 mode='lines', name='isotherm',
                                 line=dict(color='orange'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisochor), y=tscale(Tisochor), mode='lines', name='isochore',
                                 line=dict(color='gold'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale([s[0], s[0]]), y=tscale([min(Tplot), max(Tplot)]),
                                 mode='lines', name='isentrope',
                                 line=dict(color='aqua'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisenth), y=tscale(Tisenth), mode='lines', name='isenthalpe',
                                 line=dict(color='lime'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisovap), y=tscale(Tisovap), mode='lines', name='isovapour',
                                 line=dict(color='grey'), visible="legendonly"))
        fig.update_xaxes(title_text="s [kJ/kg.K]")
        fig.update_yaxes(title_text="T [°C]")
        return fig

    @output
    @render_widget
    @reactive.event(results, ignore_none=True)
    def hs() -> object:
        resultlist = results()
        try:
            s = [to_float(resultlist[var['s']])]
            h = [to_float(resultlist[var['h']])]
        except ValueError:
            s = []
            h = []
        sc = [float(resultlist[var['sc']])]
        hc = [float(resultlist[var['hc']])]

        _, _, _, _, hliq, hgas, sliq, sgas = saturation_curve()
        sisobar, _, hisobar = isobar_curve()
        _, _, sisotherm, hisotherm = isotherm_curve()
        _, _, sisochor, hisochor = isochore_curve()
        _, _, _, sisovap, hisovap = isovapour_curve()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=hscale(sliq), y=hscale(hliq),
                                 mode='lines', name='sat-liquid', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=hscale(sgas), y=hscale(hgas),
                                 mode='lines', name='sat-gas', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=hscale(sc), y=hscale(hc), mode='markers', name='critical point',
                                 marker=dict(size=10, color='black')))
        fig.add_trace(go.Scatter(x=hscale(s), y=hscale(h), mode='markers', name='result',
                                 marker=dict(size=15, color='magenta', symbol='star')))
        fig.add_trace(go.Scatter(x=hscale(sisobar), y=hscale(hisobar), mode='lines', name='isobar',
                                 line=dict(color='green'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisotherm), y=hscale(hisotherm), mode='lines', name='isotherm',
                                 line=dict(color='orange'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisochor), y=hscale(hisochor), mode='lines', name='isochore',
                                 line=dict(color='gold'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale([s[0], s[0]]), y=hscale([min(min(hliq), h[0]), max(max(hgas), h[0])]),
                                 mode='lines', name='isentrope', line=dict(color='aqua'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale([min(min(sliq), s[0]), max(max(sgas), s[0])]), y=hscale([h[0], h[0]]),
                                 mode='lines', name='isenthalpe', line=dict(color='lime'), visible="legendonly"))
        fig.add_trace(go.Scatter(x=hscale(sisovap), y=hscale(hisovap), mode='lines', name='isovapour',
                                 line=dict(color='grey'), visible="legendonly"))
        fig.update_xaxes(title_text="s [kJ/kg.K]")
        fig.update_yaxes(title_text="h [kJ/kg]")
        return fig

    @reactive.Effect()
    def _():
        method = input.method()
        if method == "MWater":
            fluid = ["H2O"]
        elif method == "PRFluid":
            fluid = PRFluidList
        elif method == "CoolProp":
            fluid = CPFluidList
        else:
            fluid = []
        ui.update_select(
            "fluid",
            choices=fluid
        )

    @reactive.Effect()
    def _():
        prop1 = input.inprop1()
        prop2 = {key: val for key, val in props.items() if key != 'u' and key != prop1}
        selected = input.inprop2()
        ui.update_select(
            "inprop2",
            choices=prop2,
            selected=selected,
        )


app = App(app_ui, server)
