from pickle import FALSE, TRUE
from networkx import DiGraph
import networkx as nx
import pandas as pd
from pathlib import Path
from math import sqrt
import numpy as np

import opendssdirect as dss

from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D

def none_to_empty(input) -> list:
    if input==None:
        out = []
    else:
        out = [input]
    return out

def term_str_to_set(term_str: str) -> set:
    term_set = set(map(int, term_str.split('.')))
    return term_set

elmt_dict = {
            'Transformer': 'xfmr', 
            'Line': 'line'
            }


class Line():
    def __init__(self, length: float, units: str, line_geo: str) -> None:
        self.length = length
        self.units = units
        self.line_geo = line_geo

    def length_in_miles(self) -> float:
        in_miles = {
            'mi': 1.,
            'kft': 1000./5280.,
            'ft': 1./5280.,
            'km': 3280.84/5280
            }
        len = self.length * in_miles[self.units]
        return len

class Load():
    def __init__(self, bus: str, kw: float, kvar: float, phases: int):
        self.kW = kw
        self.kvar = kvar
        self.phases = phases
        self.kV = 0.0
        self.terminals = ''
        self.is_delta = False

    def set_kV(self, bus_base_Vln: float):
        if self.phases==3:
            self.kV = sqrt(3)*bus_base_Vln
        else:
            self.kV = bus_base_Vln

class Xfmr():
    def __init__(self, xfmr: pd.Series) -> None:
        self.phases = int(xfmr.phases)
        self.windings = int(xfmr.windings)
        self.xpu = xfmr.X
        self.nll = xfmr.nll
        self.bus = []
        self.wnd_terms = []
        self.wnd_terms_set = []
        self.kV = []
        self.kVA = []
        self.conn = []
        self.rpu = []
        for w in range(1,(xfmr.windings+1)):
            self.bus.append(xfmr['bus'+str(w)])
            self.wnd_terms.append(str(xfmr['wnd'+str(w)+'_terms']))
            self.kV.append(xfmr['kV'+str(w)])
            self.kVA.append(xfmr['kVA'+str(w)])
            self.conn.append(xfmr['conn'+str(w)])
            self.rpu.append(float(xfmr['R']/xfmr['windings']))

    def active_terminals(self, highside: bool = False) -> set:
        if highside:
            s = 0
        else:
            s=-1
        buses = np.array(self.bus)
        windings = np.array(self.wnd_terms)
        terminals = [term_str_to_set(t) for t in windings[buses==self.bus[s]]]
        terminal_set = set().union(*terminals)
        return terminal_set

class Meter():

    def __init__(self, circ_element: str, bus1, bus2, registers: pd.DataFrame = pd.DataFrame()):
        self.name = 'meter_'+bus1+'-'+bus2
        self.c_elmt = circ_element
        self.b1 = bus1
        self.b2 = bus2
        self.registers = registers

    def to_DSS(self):
        dss_str = 'new energymeter.'+self.name
        dss_str += ' '+self.c_elmt+'.'
        dss_str += elmt_dict[self.c_elmt]+'_'+self.b1+'-'+self.b2+'_'+str(0)


class DistNetwork(DiGraph):

    def add_circuit(self, circuit_csv):
        circuit_df = pd.read_csv(circuit_csv)
        circuit_dict = circuit_df.squeeze().to_dict()
        self.circuit = circuit_dict    
    
    def add_nodes(self, node_csv):
        nodes_df = pd.read_csv(node_csv)
        for bus in nodes_df.index:
            bus_name = nodes_df['Bus'][bus]
            Loc_x = nodes_df['Graph_Loc_x'][bus]
            Loc_y = nodes_df['Graph_Loc_y'][bus]
            graph_coords = (Loc_x, Loc_y)
            self.add_node(bus_name, coords = graph_coords, hot_terminals = set())
        self.nodes[self.circuit['slack_bus']]['hot_terminals'].update({1,2,3})
    
    def add_line_geoms(self, lineGeoms_csv):
        line_geoms_df = pd.read_csv(lineGeoms_csv, 
                                    index_col='LineGeo'
                                    )
        self.line_geoms = line_geoms_df
        phase_dict = {'a':'1',
                      'b':'2',
                      'c':'3',
                      'ab': '1.2',
                      'ac': '1.3',
                      'bc': '2.3',
                      'abc': '1.2.3'}
        self.line_geoms['cond_pos'] = self.line_geoms['cond_pos'].map(phase_dict)
        self.line_geoms['terminals'] = self.line_geoms['cond_pos'].map(term_str_to_set)


    def get_network_wires(self) -> set:
        OH_line_geos = self.line_geoms['Type']=='OH'
        n_cond_select = self.line_geoms['Nphases'] != self.line_geoms['Nconds']
        ph_conductors = set(self.line_geoms['phase_cond'][OH_line_geos])
        n_conductors = set(self.line_geoms['neutral_cond'][n_cond_select])
        return n_conductors.union(ph_conductors)

    def get_network_cnwires(self) -> set:
        CN_line_geos = self.line_geoms['Type']=='CN'
        cn_conductors = set(self.line_geoms['neutral_cond'][CN_line_geos])
        return cn_conductors
    
    def add_wiredata(self, wireData_csv):
        columns = ['Type',
                   'GMRac',
                   'GMRunits',
                   'rac',
                   'runits',
                   'normamps'
                   ]
        wire_data_df = pd.read_csv(wireData_csv, index_col='Type', usecols=columns)
        wires = wire_data_df.index
        wset = wires.isin(self.get_network_wires())
        self.wire_data = wire_data_df.loc[wset]

    def add_CNdata(self, CNData_csv):
        columns = ['Type',
                   'GMRac',
                   'GmrStrand',
                   'GMRunits',
                   'Rac',
                   'Rstrand',
                   'Runits',
                   'DiaCable',
                   'DiaIns',
                   'DiaStrand',
                   'diam',
                   'Radunits',
                   'normamps',
                   'k'
                   ]
        CN_data_df = pd.read_csv(CNData_csv, index_col='Type', usecols=columns)
        cn_wires = CN_data_df.index
        cnset = cn_wires.isin(self.get_network_cnwires())
        self.CN_data = CN_data_df.loc[cnset]
    
    def add_edge_DN(self, u, v, 
                    line: Line = None, 
                    xfmr: Xfmr = None, 
                    length: float = 0.):
       l = none_to_empty(line)
       x = none_to_empty(xfmr)
       self.add_edge(u, v, lines = l, xfmrs = x, length = length)

    def add_lines(self, lines_csv):
        lines_df = pd.read_csv(lines_csv)
        for line in lines_df.index:
            Bus1 = lines_df['bus1'][line]
            Bus2 = lines_df['bus2'][line]
            length = lines_df['Length'][line]
            units = lines_df['Units'][line]
            line_geo = lines_df['LineGeometry'][line]
            newline = Line(length, units, line_geo)
            hot_terminals = set(self.line_geoms['terminals'][line_geo])
            self.nodes[Bus2]['hot_terminals'].update(hot_terminals.difference({0}))
            if self.has_edge(Bus1, Bus2):
                self[Bus1][Bus2]['lines'].append(newline)
                self[Bus1][Bus2]['length'] = newline.length_in_miles()
                if self[Bus1][Bus2]['xfmrs']:
                    print(f'Adding line across pre-existing xfmr betwen {Bus1} and {Bus2}')
                else:
                    pass
            else:
                self.add_edge_DN(
                          Bus1, 
                          Bus2,
                          line = newline,
                          length = newline.length_in_miles(),
                          )
    
    def add_loads(self, loads_csv):
        loads_df = pd.read_csv(loads_csv)
        for load in loads_df.index:
            bus = loads_df['Bus'][load]
            load_kW = loads_df['kW'][load]
            load_kVAr = loads_df['kVAr'][load]
            phases = loads_df['Phases'][load]
            newload = Load(bus=bus, kw=load_kW, kvar=load_kVAr, phases=phases)

            ht = self.nodes[bus]['hot_terminals']
            load_terms = str(loads_df['CondPos'][load])
            load_term_set = term_str_to_set(load_terms)
            for t in load_term_set:
                if t not in ht:
                    print(f'Load at {bus} connected to terminal {t} which is not energized')
                else:
                    pass
            newload.terminals = load_terms
            bus_base_Vln = self.nodes[bus]['Vln_base']
            newload.set_kV(bus_base_Vln)
            if newload.phases>1:
                upstrm_xfmr = self.get_upstream_xfmr(bus)
                if upstrm_xfmr.conn[-1]=='delta':
                    newload.is_delta = True
            self.nodes[bus]['load'] = newload

    def downstream_loads(self, node: str, as_kVA: bool = False):
        ds_nodes = [node]
        ds_nodes += list(nx.descendants(self, node))
        ds_kW = 0.0
        ds_kVAr = 0.0
        for n in ds_nodes:
            if 'load' in self.nodes[n]:
                ds_kW += self.nodes[n]['load'].kW
                ds_kVAr += self.nodes[n]['load'].kvar
            else:
                pass
        if as_kVA:
            kVA = sqrt(ds_kW**2 + ds_kVAr**2)
            return kVA
        else:
            load_tuple = (ds_kW, ds_kVAr)
            return load_tuple

    def add_xfmrs(self, xfmr_csv):
        xfmr_df = pd.read_csv(xfmr_csv)
        for x in xfmr_df.index:
            xfmr = xfmr_df.loc[x]
            new_xfmr = Xfmr(xfmr)
            Bus1 = new_xfmr.bus[0]
            Bus2 = new_xfmr.bus[-1]
            self.nodes[Bus2]['hot_terminals'].update(new_xfmr.active_terminals())
            if self.has_edge(Bus1, Bus2):
                self[Bus1][Bus2]['xfmr'].append(new_xfmr)
                self[Bus1][Bus2]['length'] = 0.
                if self[Bus1][Bus2]['lines']:
                    print (f'Adding xfmr across pre-existing line betwen {Bus1} and {Bus2}')
                else:
                    pass
            else:
                self.add_edge_DN(Bus1, Bus2, xfmr = new_xfmr)

    def calc_electrical_distance(self):
        source_node = self.circuit['slack_bus']
        for n in self.nodes:
            len = nx.dijkstra_path_length(self, source=source_node, target=n , weight='length')
            self.nodes[n]['dist_from_source'] = len

    def get_upstream_xfmr(self, node: str) -> Xfmr:
        b1 = node
        while b1 != self.circuit['slack_bus']:
            b2 = next(self.predecessors(b1))
            if self[b2][b1]['xfmrs']: #If xfmrs list is non-empty, then:
                xfmr = self[b2][b1]['xfmrs'][0]
                return xfmr
            else: # if xfmrs list is empty then:
                b1 = b2
    
    def calc_Vbases(self):
        for n in self.nodes:
            if n==self.circuit['slack_bus']:
                Vbase = (1/sqrt(3)) * self.circuit['basekv']
            else:
                upstrm_xfmr = self.get_upstream_xfmr(n)
                if upstrm_xfmr.phases==3:
                    Vbase = (1/sqrt(3)) * upstrm_xfmr.kV[-1]
                elif upstrm_xfmr.conn[-1]=='delta':
                    Vbase = (1/sqrt(3)) * upstrm_xfmr.kV[-1]
                else:
                    Vbase = upstrm_xfmr.kV[-1]
            self.nodes[n]['Vln_base'] = Vbase*1000

    def get_line_terminals(self, node: str):
        prev_node = next(self.predecessors(node))

    def check_floating_lines(self):
        for n in self.nodes:
            n_ht = self.nodes[n]['hot_terminals']
            for m in self.predecessors(n):
                m_ht = self.nodes[m]['hot_terminals']
                for t in n_ht:
                    if t not in m_ht:
                        print(f'There exists no terminal {t} connection between nodes {n} and {m}')
                    else:
                        pass

    def add_network_components(self, component_csvs: dict):
        
        circuit_csv = Path('..', 'network-data', component_csvs['circuit'])
        self.add_circuit(circuit_csv)
        
        bus_csv = Path('..', 'network-data', component_csvs['nodes'])
        self.add_nodes(bus_csv)
        
        lineGeo_csv = Path('..', 'network-data', component_csvs['line_geos'])
        self.add_line_geoms(lineGeo_csv)
        
        wireData_csv = Path('..', 'network-data', component_csvs['wire_data'])
        self.add_wiredata(wireData_csv)
        CNData_csv = Path('..', 'network-data', component_csvs['CN_data'])
        self.add_CNdata(CNData_csv)
        
        lines_csv = Path('..', 'network-data', component_csvs['lines'])
        self.add_lines(lines_csv)

        xfmrs_csv = Path('..', 'network-data', component_csvs['xfmrs'])
        self.add_xfmrs(xfmrs_csv)

        self.calc_Vbases()

        loads_csv = Path('..', 'network-data', component_csvs['loads'])
        self.add_loads(loads_csv)

        self.calc_electrical_distance()
        
        self.check_floating_lines()



    def clear_DSS(self):
        dss.run_command("clear")

    def solve_DSS(self):
        dss.run_command('set maxiterations=300')
        dss.run_command('set maxcontroliter=300')
        dss.run_command('solve')

    def new_circuit_DSS(self):
        cmd = dss.run_command
        circ = self.circuit

        dss_str = 'new circuit.slack bus1=' + circ['slack_bus']
        dss_str += ' basekv=' + str(circ['basekv'])
        dss_str += ' pu=' + str(circ['pu'])
        dss_str += ' phases=' + str(circ['phases'])
        dss_str += ' Z1=['+str(circ['R_1'])+', '+str(circ['X_1'])+']' #Z1=[R_1, X_1]
        dss_str += ' Z0=['+str(circ['R_0'])+', '+str(circ['X_0'])+']' #Z0=[R_0, X_0]
        cmd(dss_str)

    def new_wire_DSS(self, wire_name: str, wire_data: pd.Series):
        dss_str = 'new wiredata.' + wire_name
        for index, value in wire_data.items():
            dss_str += ' ' + index + '=' + str(value)
        dss.run_command(dss_str)

    def new_cnwire_DSS(self, cnwire_name: str, cnwire_data: pd.Series):
        dss_str = 'new cndata.' + cnwire_name
        for index, value in cnwire_data.items():
            dss_str += ' ' + index + '=' + str(value)
        dss.run_command(dss_str)

    def new_linegeo_DSS(self, lg_name: str, lg: pd.Series):
        dss_str = 'new lineGeometry.' + lg_name
        dss_str += ' Nconds=' + str(lg.Nconds)
        dss_str += ' Nphases=' + str(lg.Nphases)
        dss_str += ' Units=' + str(lg.Units)
        for cond in range(1, lg.Nconds+1):
            dss_str += ' Cond=' + str(cond)
            if lg.Type=='OH':
                dss_str += ' wire=' + str(lg.phase_cond)
            elif lg.Type=='CN':
                dss_str += ' cnCable=' + str(lg.neutral_cond)
            else:
                pass
            xcol = str(cond) + '_x'
            zcol = str(cond) + '_z'
            dss_str += ' X=' + str(lg[xcol])
            dss_str += ' H=' + str(lg[zcol])
            if lg.Nconds > lg.Nphases:
                dss_str += ' reduce=y'
            else:
                dss_str += ' reduce=n'
        dss.run_command(dss_str)

    def new_line_DSS(self, b1: str, b2: str, line: Line, id: int):
        dss_str = 'New Line.line_'+b1+'-'+b2+'_'+str(id)
        dss_str += ' Length=' + str(line.length)
        dss_str += ' Units=' + line.units
        lg = self.line_geoms.loc[line.line_geo]
        dss_str += ' Bus1=' + b1 + '.' + str(lg.cond_pos)
        dss_str += ' Bus2=' + b2 + '.' + str(lg.cond_pos)
        dss_str += ' phases=' + str(lg.Nphases)
        dss_str += ' Geometry=' + str(line.line_geo)
        dss.run_command(dss_str)

    def new_xfmr_DSS(self, b1: str, b2: str, xfmr: Xfmr, id: int):
        dss_str = 'new transformer.xfmr_'+b1+'-'+b2+'_'+str(id)
        dss_str += ' phases=' + str(xfmr.phases)
        dss_str += ' windings=' + str(xfmr.windings)
        dss_str += ' XHL=' + str(xfmr.xpu)
        dss_str += ' %noloadloss=' + str(xfmr.nll)
        for w in range(xfmr.windings):
            dss_str += ' wdg=' + str(w+1)
            dss_str += ' bus=' + xfmr.bus[w] + '.' + str(xfmr.wnd_terms[w])
            dss_str += ' conn=' + str(xfmr.conn[w])
            dss_str += ' kV=' + str(xfmr.kV[w])
            dss_str += ' kVA=' + str(xfmr.kVA[w])
            dss_str += ' %R=' + str(xfmr.rpu[w])
        dss.run_command(dss_str)


    def new_load_DSS(self, node: str, load: Load):
        dss_str = 'New Load.load_' + node
        dss_str += ' Bus1=' + node + '.' + load.terminals
        dss_str += ' Phases=' + str(load.phases)
        dss_str += ' kw=' + str(load.kW)
        dss_str += ' kvar=' + str(load.kvar)
        dss_str += ' kv=' + str(load.kV)
        if load.is_delta:
            dss_str += ' conn=delta'
        dss.run_command(dss_str)


    def compile_DSS(self):
        self.clear_DSS()

        self.new_circuit_DSS()

        #compile wriedata
        for wire in self.wire_data.index:
            wd = self.wire_data.loc[wire]
            self.new_wire_DSS(wire, wd)
        
        #compile CNdata
        for cnwire in self.CN_data.index:
            cnwd = self.CN_data.loc[cnwire]
            self.new_cnwire_DSS(cnwire, cnwd)

        #compile line geometries
        for line_geo in self.line_geoms.index:
            lg = self.line_geoms.loc[line_geo]
            self.new_linegeo_DSS(line_geo, lg)
        
        #compile lines
        for b1, b2, e in self.edges.data('lines'):
            if e: #in lines attr is a non-empty list then:
                for l in range(len(self[b1][b2]['lines'])):
                    line = self[b1][b2]['lines'][l]
                    self.new_line_DSS(b1, b2, line, l)
            else: #if lines attr is an empty list then:
                pass

        #compile xfmrs
        xfmr_Vbases = set()
        for b1, b2, e in self.edges.data('xfmrs'):
            if e: #if xfmrs attr is a non-empty list then:
                for x in range(len(self[b1][b2]['xfmrs'])):
                    xfmr = self[b1][b2]['xfmrs'][x]
                    xfmr_Vbases.union(set(xfmr.kV))
                    self.new_xfmr_DSS(b1, b2, xfmr=xfmr, id=x)
            else: #if xfmrs attr is an empty list then:
                pass
        vb_str = 'Set Voltagebases = '
        vb_str += '[' + ', '.join([str(v) for v in xfmr_Vbases]) + ']'
        dss.run_command(vb_str)
        dss.run_command('Calcvoltagebases')

        #compile loads
        for n, load in self.nodes(data='load'):
            if type(load)==Load:
                self.new_load_DSS(n, load)
            else:
                pass

    def record_DSS_bus_voltages(self):
        for bus in self.nodes:
            dss.Circuit.SetActiveBus(bus)
            voltage = np.array(dss.Bus.VMagAngle())
            vpu = np.array(dss.Bus.PuVoltage())
            voltage = voltage.reshape(int(voltage.size/2),2)
            vpu = vpu.reshape(int(vpu.size/2),2)
            vpu = np.array([n[0] + n[1]*1j for n in vpu])
            terms = self.nodes[bus]['hot_terminals']
            vbase = self.nodes[bus]['Vln_base']
            bus_voltage = pd.DataFrame({
                                        'Vmag': voltage[:,0], 
                                        'Vang': voltage[:, 1], 
                                        'Vpu': (1/(vbase))*voltage[:,0],
                                        'Vpu_dss': abs(vpu)
                                        }, index=terms)
            self.nodes[bus]['voltage'] = bus_voltage

    def plot_bus_voltages(self, source = False):
        if source:
            src = source
        else:
            src = self.circuit['slack_bus']
        subgraph = self.subgraph(nx.descendants(self, src).union({src}))
        pltlines = np.empty((1,2,2))
        clist = []
        cdict = {1: 'r', 2: 'g', 3: 'b'}
        for b1, b2 in subgraph.edges:
            for t in self.nodes[b2]['hot_terminals']:
                pltline = np.empty((1,2,2))
                pltline[0, 0, 0] = self.nodes[b1]['dist_from_source']
                pltline[0, 1, 0] = self.nodes[b2]['dist_from_source']
                pltline[0, 0, 1] = self.nodes[b1]['voltage'].Vpu[t]
                pltline[0, 1, 1] = self.nodes[b2]['voltage'].Vpu[t]
                pltlines = np.append(pltlines, pltline, axis=0)
                clist.append(cdict[t])
        
        branchtext = []
        for n in subgraph.nodes:
            if len(DN[n])==0:
                for t in self.nodes[n]['hot_terminals']:
                    dist = self.nodes[n]['dist_from_source']
                    vpu = self.nodes[n]['voltage'].Vpu[t]
                    branchtext.append((dist, vpu, n+'.'+str(t)))

        fig, ax = plt.subplots(figsize=(12,10))
        dmin = pltlines[1:, :, 0].min()
        dmax = pltlines[1:, :, 0].max()
        vmin = pltlines[1:, :, 1].min()
        vmax = pltlines[1:, :, 1].max()
        vdelta = vmax - vmin
        ax.set_xlim(dmin,dmax+0.1*(dmax-dmin))
        ax.set_ylim(vmin - 0.1*vdelta, vmax + 0.1*vdelta)

        ax.set_xlabel('distance from substation (mi)')
        ax.set_ylabel('V_ln (pu)')

        cdict_inv = {v: k for k, v in cdict.items()}
        phase_dict = {1: 'phase a', 2: 'phase b', 3: 'phase c'}
        leg_lines = []
        leg_labels = []
        for t in sorted([cdict_inv[c] for c in set(clist)]):
            leg_lines.append(Line2D([0], [0], color=cdict[t], lw=3))
            leg_labels.append(phase_dict[t])
        ax.legend(leg_lines, leg_labels, loc='lower left')
            

        for t in branchtext:
            plt.text(t[0], t[1], t[2])

        line_segments = LineCollection(pltlines[1:, :, :], colors = clist)
        ax.add_collection(line_segments)

        plt.show()