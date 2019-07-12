'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from . import order
from .order import Order
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from matplotlib.backends.backend_pdf import PdfPages

LOGGER = log.LOGGER

class Density(SysInfo):

    def __init__(self,inputfilename="inputfile"):
        super().__init__(inputfilename)
    
        self.lipid_type_items = ' '.join(self.molecules)
        self.u = mda.Universe(self.gropath,self.trjpath)
        
        self.main_lipid = ''.join(self.molecules[0])
        self.sterol_lipid = ''.join(self.molecules[1])
        
        self.lipid_type_items = ' '.join(self.molecules)
        
        self.lipid_head = self.u.select_atoms('resname {} and name P'.format(self.main_lipid)) 
        
        head_atoms = lipidmolecules.head_atoms_of(self.sterol_lipid)
        head_atom = head_atoms[0]

        self.sterol_head = self.u.select_atoms('resname {} and name {}'.format(self.sterol_lipid,head_atom)) 
        
    def density_main_lipid(self,start_time,end_time):

        cwd = os.getcwd()
        index_density_main_lipid = ''.join(cwd + '/' + 'index_density_main_lipid.ndx')

        density_main_lipid_raw = ''.join(cwd + '/' + 'density_main_lipid_raw')

        # get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_density, \
        #     '-select', '(resname {} and name P)'.format(self.main_lipid)]
        # out, err = exec_gromacs(get_selection)

        with mda.selections.gromacs.SelectionWriter('index_density_main_lipid.ndx', mode='w') as ndx:
            ndx.write(self.lipid_head, name='main_lipid')

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_main_lipid, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_main_lipid_raw, '-center', '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        with open("density_main.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("density_main_lipid_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_main_lipid.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))


    def density_sterol_lipid(self,start_time,end_time):

        cwd = os.getcwd()
        index_density_sterol_lipid = ''.join(cwd + '/' + 'index_density_sterol_lipid.ndx')

        density_sterol_lipid_raw = ''.join(cwd + '/' + 'density_sterol_lipid_raw')

        # get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_density, \
        #     '-select', '(resname {} and name P)'.format(self.main_lipid)]
        # out, err = exec_gromacs(get_selection)

        with mda.selections.gromacs.SelectionWriter('index_density_sterol_lipid.ndx', mode='w') as ndx:
            ndx.write(self.sterol_head, name='sterol_lipid')

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_sterol_lipid, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_sterol_lipid_raw, '-center', '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        with open("density_sterol.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("density_sterol_lipid_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_sterol_lipid.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
