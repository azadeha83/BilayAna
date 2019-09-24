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

    def rdf(self,ref, sel, start_time, end_time, step, seltype='atom', selrpos='atom', binsize=0.002):
        
        n_ref = ref.split(' ')[1]
        n_sel = sel.split(' ')[1]
        
        file_name = n_ref + '_' + n_sel

        ref_sel = self.u.select_atoms('{}'.format(ref))
        sel_sel = self.u.select_atoms('{}'.format(sel))
        
        ref_sel_center = ref_sel.center_of_geometry()        
        sel_sel_center = sel_sel.center_of_geometry()        
        
        selection1 = self.u.select_atoms('{} and (prop z > {})'.format(ref,ref_sel_center[2])).resids
        selection2 = self.u.select_atoms('{} and (prop z > {})'.format(sel,sel_sel_center[2])).resids
        
        cwd = os.getcwd()
        rdf_raw = ''.join(cwd + '/' + 'rdf_raw_' + '{}'.format(file_name))
        rdf_cum_raw = ''.join(cwd + '/' + 'rdf_cum_raw_' + '{}'.format(file_name))
        
        '''
        rdf_index = ''.join(cwd + '/' + 'rdf_index_' + '{}'.format(file_name) + '.ndx')

        with mda.selections.gromacs.SelectionWriter('rdf_index_' + '{}'.format(file_name) + '.ndx', mode='w') as ndx:
            ndx.write(selection1, name='{}'.format(ref))
            ndx.write(selection2, name='{}'.format(sel))
        '''
        get_rdf = [GMXNAME, 'rdf', '-f', self.trjpath, '-s', self.tprpath, \
            '-o', rdf_raw, '-cn', rdf_cum_raw, '-ref', '{}'.format(ref)+' and resid '+' '.join(map(str, selection1)), \
                '-sel', '{}'.format(sel)+' and resid '+' '.join(map(str, selection2)), '-xy',
                '-selrpos', selrpos, '-seltype', seltype, '-bin', str(binsize), '-b', str(start_time), '-e', str(end_time), '-xvg', 'none','-dt', str(step)]
                        
        out, err = exec_gromacs(get_rdf)
        
        with open("rdf.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

    
    
