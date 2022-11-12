from math import sqrt
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from networkx import DiGraph


class Line():
    def __init__(self, length: float, units: str, line_geo: str) -> None:
        self.length = length
        self.units = units
        self.line_geo = line_geo

class Load():
    def __init__(self, kw: float, kvar: float):
        self.kW = kw
        self.kvar = kvar

class Xfmr():
    def __init__(self, 
                kV_hv: float, 
                kV_lv: float, 
                kVA_hv: float, 
                kVA_lv: float, 
                conn_hv: str, 
                conn_lv: str,
                phases: int,
                rpu: float,
                xpu: float
                ) -> None:
        self.kV_hv = kV_hv
        self.kV_lv = kV_lv
        self.kVa_hv = kVA_hv
        self.kVA_lv = kVA_lv
        self.conn_hv = conn_hv
        self.conn_lv = conn_lv
        self.phases = phases
        self.rpu = rpu
        self.xpu = xpu

class DistNetwork(DiGraph):

    def Add_Nodes(self, node_csv):
        nodes_df = pd.read_csv(node_csv)
        for bus in nodes_df.index:
            bus_name = nodes_df['Bus'][bus]
            Loc_x = nodes_df['Graph_Loc_x'][bus]
            Loc_y = nodes_df['Graph_Loc_y'][bus]
            graph_coords = (Loc_x, Loc_y)
            self.add_node(bus_name, coords = graph_coords)
    
    def Add_LineGeoms(self, lineGeoms_csv):
        line_geoms_df = pd.read_csv(lineGeoms_csv, 
                                    index_col='LineGeo'
                                    )
        self.line_geoms = line_geoms_df
        phase_dict = {'a':'.1',
                      'b':'.2',
                      'c':'.3',
                      'ab': '.1.2',
                      'ac': '.1.3',
                      'bc': '.2.3',
                      'abc': '.1.2.3'}
        self.line_geoms['cond_pos'] = self.line_geoms['cond_pos'].map(phase_dict)


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
    
    def Add_WireData(self, wireData_csv):
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

    def Add_CNData(self, CNData_csv):
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

    def Add_Lines(self, lines_csv):
        lines_df = pd.read_csv(lines_csv)
        for line in lines_df.index:
            Bus1 = lines_df['bus1'][line]
            Bus2 = lines_df['bus2'][line]
            length = lines_df['Length'][line]
            units = lines_df['Units'][line]
            line_geo = lines_df['LineGeometry'][line]
            newline = Line(length, units, line_geo)
            if self.has_edge(Bus1, Bus2):
                self[Bus1][Bus2]['line'] = newline
            else:
                self.add_edge(
                          Bus1, 
                          Bus2,
                          line = newline,
                          )
    
    def Add_Loads(self, loads_csv):
        loads_df = pd.read_csv(loads_csv)
        for load in loads_df.index:
            bus = loads_df['Node'][load]
            load_kW = loads_df['kW'][load]
            load_kVAr = loads_df['kVAr'][load]
            newload = Load(load_kW, load_kVAr)
            self.nodes[bus]['load'] = newload

    def downstream_loads(self, node: str, as_kVA: bool = False):
        ds_nodes = [node]
        ds_nodes += list(nx.descendants(self, node))
        ds_kW = 0.0
        ds_kVAr = 0.0
        for n in ds_nodes:
            if 'load_node' in self.nodes[n]:
                ds_kW += self.nodes[n]['load_kW']
                ds_kVAr += self.nodes[n]['load_kVAr']
            else:
                pass
        if as_kVA:
            kVA = sqrt(ds_kW**2 + ds_kVAr**2)
            return kVA
        else:
            load_tuple = (ds_kW, ds_kVAr)
            return load_tuple

    def Add_Xfmrs(self, xfmr_csv):
        xfmr_df = pd.read_csv(xfmr_csv)
        #print(xfmr_df)
        for xfmr in xfmr_df.index:
            Bus1 = xfmr_df['bus1_hv'][xfmr]
            Bus2 = xfmr_df['bus2_lv'][xfmr]
            xfmr_details = {
                'kV_hv': xfmr_df['kV_hv'][xfmr],
                'kV_lv': xfmr_df['kV_lv'][xfmr],
                'kVA_hv': xfmr_df['kVA_hv'][xfmr],
                'kVA_lv': xfmr_df['kVA_hv'][xfmr],
                'conn_hv': xfmr_df['conn_hv'][xfmr],
                'conn_lv': xfmr_df['conn_lv'][xfmr],
                'phases': xfmr_df['phases'][xfmr],
                'rpu': xfmr_df['R'][xfmr],
                'xpu': xfmr_df['X'][xfmr]
            }
            new_xfmr = Xfmr(*xfmr_details)
            if self.has_edge(Bus1, Bus2):
                self[Bus1][Bus2]['xfmr'] = new_xfmr
            else:
                self.add_edge(Bus1, Bus2, xfmr = new_xfmr)

    def add_network_components(self, component_csvs: dict):
        
        bus_csv = Path('..', 'network-data', component_csvs['nodes'])
        self.Add_Nodes(bus_csv)
        
        lineGeo_csv = Path('..', 'network-data', component_csvs['line_geos'])
        self.Add_LineGeoms(lineGeo_csv)
        
        wireData_csv = Path('..', 'network-data', component_csvs['wire_data'])
        self.Add_WireData(wireData_csv)
        CNData_csv = Path('..', 'network-data', component_csvs['CN_data'])
        self.Add_CNData(CNData_csv)
        
        lines_csv = Path('..', 'network-data', component_csvs['lines'])
        self.Add_Lines(lines_csv)
        
        loads_csv = Path('..', 'network-data', component_csvs['loads'])
        self.Add_Loads(loads_csv)

        xfmrs_csv = Path('..', 'network-data', component_csvs['xfmrs'])
        self.Add_Xfmrs(xfmrs_csv)