'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

LOGGER = log.LOGGER

class Rdf(SysInfo):

    def __init__(self,inputfilename="inputfile"):
        super().__init__(inputfilename)
    
        self.u = mda.Universe(self.gropath,self.trjpath)

    def rdf(self, ref, sel, ref0, sel0, start_time, end_time, dt, seltype='atom', selrpos='atom', binsize=0.002):
        
        w = mda.Universe(self.gropath)
        
        file_name = ref.split(' ')[1] + '_' + sel.split(' ')[1]

        ref_sel = w.select_atoms('{}'.format(ref))
        sel_sel = w.select_atoms('{}'.format(sel))
        
        ref_sel_center = ref_sel.center_of_geometry()        
        sel_sel_center = sel_sel.center_of_geometry()        
        
        selection1 = list(set(w.select_atoms('{} and (prop z > {})'.format(ref0,ref_sel_center[2])).resids))
        selection2 = list(set(w.select_atoms('{} and (prop z > {})'.format(sel0,sel_sel_center[2])).resids))
        
        print(selection1)
        print(selection2)
        cwd = os.getcwd()
        rdf_raw = ''.join(cwd + '/' + 'rdf_raw_' + '{}'.format(file_name)+ '_' +'{}'.format(ref.split(' ')[-1])+ '_' +'{}'.format(sel.split(' ')[-1]))
        rdf_cum_raw = ''.join(cwd + '/' + 'rdf_cum_raw_' + '{}'.format(file_name)+ '_' +'{}'.format(ref.split(' ')[-1])+ '_' +'{}'.format(sel.split(' ')[-1]))
        
        '''
        rdf_index = ''.join(cwd + '/' + 'rdf_index_' + '{}'.format(file_name) + '.ndx')

        with mda.selections.gromacs.SelectionWriter('rdf_index_' + '{}'.format(file_name) + '.ndx', mode='w') as ndx:
            ndx.write(selection1, name='{}'.format(ref))
            ndx.write(selection2, name='{}'.format(sel))
        '''
        get_rdf = [GMXNAME, 'rdf', '-f', self.trjpath, '-s', self.tprpath, \
            '-o', rdf_raw, '-cn', rdf_cum_raw, '-ref', '{}'.format(ref)+' and resid '+' '.join(map(str, selection1)), \
                '-sel', '{}'.format(sel)+' and resid '+' '.join(map(str, selection2)), '-xy',
                '-selrpos', selrpos, '-seltype', seltype, '-bin', str(binsize), '-b', str(start_time), '-e', str(end_time), '-dt', str(dt), '-xvg', 'none']
                        
        out, err = exec_gromacs(get_rdf)
        
        with open("rdf.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

    
    
