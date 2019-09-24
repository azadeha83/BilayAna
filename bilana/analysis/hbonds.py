'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import os
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
import MDAnalysis as mda
import pandas as pd
import numpy as np
import csv

LOGGER = log.LOGGER

class Hbonds(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)

            self.lipid_types_mainlipid = ''.join(self.molecules[0])
            self.lipid_types_sterol = ''.join(self.molecules[1])

    def hbonds_sterol_lipid(self):
        ''' This module calculates the hbonds between sterol and lipid molecules '''

        cwd = os.getcwd()
       
        hbond_raw = ''.join(cwd + '/' + 'hbond_sterol_lipid_raw')
        hbond_log = ''.join(cwd + '/' + 'hbond_sterol_lipid_log')
        index_hbond = ''.join(cwd + '/' + 'index_hbond.ndx')
        
        get_selection = [GMXNAME, 'select', '-f', self.gropath, '-s', self.tprpath, '-on', index_hbond, \
            '-select', '(resname {} {})'.format(self.lipid_types_mainlipid,self.lipid_types_sterol)]
        
        out, err = exec_gromacs(get_selection)
        
        get_hbond = [GMXNAME, 'hbond', '-f', self.trjpath, '-s', self.tprpath, '-n', index_hbond, '-ac', hbond_raw, \
            '-g', hbond_log, '-num', hbond_raw]
        
        out, err = exec_gromacs(get_hbond)

        with open("hbond_sterol_lipid_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("hbond_sterol_lipid.dat" , 'w') as fout:

            fout.write('Time\tn_hbonds\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))

    def hbonds_sterol_sterol(self):
        ''' This module calculates the hbonds between sterol and sterol molecules '''

        cwd = os.getcwd()
        
        hbond_raw = ''.join(cwd + '/' + 'hbond_sterol_sterol_raw')
        hbond_log = ''.join(cwd + '/' + 'hbond_sterol_sterol_log')
        index_hbond = ''.join(cwd + '/' + 'index_hbond_sterol.ndx')
        
        get_selection = [GMXNAME, 'select', '-f', self.gropath, '-s', self.tprpath, '-on', index_hbond, \
            '-select', '(resname {})'.format(self.lipid_types_sterol)]
        
        out, err = exec_gromacs(get_selection)

        get_hbond = [GMXNAME, 'hbond', '-f', self.trjpath, '-s', self.tprpath, '-ac', hbond_raw, \
            '-g', hbond_log, '-num', hbond_raw]
                   
        out, err = exec_gromacs(get_hbond)

        with open("hbond_sterol_sterol_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("hbond_sterol_sterol.dat" , 'w') as fout:

            fout.write('Time\tn_hbonds\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))