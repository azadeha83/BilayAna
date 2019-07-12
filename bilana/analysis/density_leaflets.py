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
        #print(len(self.u.trajectory))

        w = mda.Universe(self.gropath)
        
        self.main_lipid = ''.join(self.molecules[0])
        self.sterol_lipid = ''.join(self.molecules[1])
        
        self.lipid_type_items = ' '.join(self.molecules)
        
        main_lipid_selection = w.select_atoms('resname {} and name P'.format(self.main_lipid))
        
        main_lipid_selection_center = main_lipid_selection.center_of_geometry()        

        lipid_head_up = self.u.select_atoms('resname {} and name P and (prop z > {})'.format(self.main_lipid,main_lipid_selection_center[2])) 
        lipid_head_down = self.u.select_atoms('resname {} and name P and (prop z < {})'.format(self.main_lipid,main_lipid_selection_center[2]))
        
        lipid_head_up_resids = lipid_head_up.resids        
        lipid_head_down_resids = lipid_head_down.resids
        
        self.lipid_head_up = self.u.select_atoms('resname {} and resid {}'.format(self.main_lipid,' '.join(str(v) for v in lipid_head_up_resids))) 
        self.lipid_head_down = self.u.select_atoms('resname {} and resid {}'.format(self.main_lipid,' '.join(str(v) for v in lipid_head_down_resids)))
        
        sterol_head_up = self.u.select_atoms('resname {} and name O3 and (prop z > {})'.format(self.sterol_lipid,main_lipid_selection_center[2])) 
        sterol_head_down = self.u.select_atoms('resname {} and name O3 and (prop z < {})'.format(self.sterol_lipid,main_lipid_selection_center[2])) 
        
        sterol_head_up_resids = sterol_head_up.resids
        sterol_head_down_resids = sterol_head_down.resids
        
        self.sterol_head_up = self.u.select_atoms('resname {} and resid {}'.format(self.sterol_lipid,' '.join(str(v) for v in sterol_head_up_resids))) 
        self.sterol_head_down = self.u.select_atoms('resname {} and resid {}'.format(self.sterol_lipid,' '.join(str(v) for v in sterol_head_down_resids)))

    def density_main_lipid(self,start_time,end_time):

        cwd = os.getcwd()
        index_density_main_lipid_upper = ''.join(cwd + '/' + 'index_density_main_lipid_upper.ndx')
        index_density_main_lipid_lower = ''.join(cwd + '/' + 'index_density_main_lipid_lower.ndx')

        density_main_lipid_upper_raw = ''.join(cwd + '/' + 'density_main_lipid_upper_raw')
        density_main_lipid_lower_raw = ''.join(cwd + '/' + 'density_main_lipid_lower_raw')

        # get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_density, \
        #     '-select', '(resname {} and name P)'.format(self.main_lipid)]
        # out, err = exec_gromacs(get_selection)

        with mda.selections.gromacs.SelectionWriter('index_density_main_lipid_upper.ndx', mode='w') as ndx:
            ndx.write(self.lipid_head_up, name='upper_main_lipid')
        
        with mda.selections.gromacs.SelectionWriter('index_density_main_lipid_lower.ndx', mode='w') as ndx:
            ndx.write(self.lipid_head_down, name='lower_main_lipid')

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_main_lipid_upper, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_main_lipid_upper_raw, '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_main_lipid_lower, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_main_lipid_lower_raw, '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        with open("density_main.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("density_main_lipid_upper_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_main_lipid_upper.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))

        with open("density_main_lipid_lower_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_main_lipid_lower.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))

    def density_sterol_lipid(self,start_time,end_time):

        cwd = os.getcwd()
        index_density_sterol_lipid_upper = ''.join(cwd + '/' + 'index_density_sterol_lipid_upper.ndx')
        index_density_sterol_lipid_lower = ''.join(cwd + '/' + 'index_density_sterol_lipid_lower.ndx')

        density_sterol_lipid_upper_raw = ''.join(cwd + '/' + 'density_sterol_lipid_upper_raw')
        density_sterol_lipid_lower_raw = ''.join(cwd + '/' + 'density_sterol_lipid_lower_raw')

        # get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_density, \
        #     '-select', '(resname {} and name P)'.format(self.main_lipid)]
        # out, err = exec_gromacs(get_selection)

        with mda.selections.gromacs.SelectionWriter('index_density_sterol_lipid_upper.ndx', mode='w') as ndx:
            ndx.write(self.sterol_head_up, name='upper_sterol_lipid')
        
        with mda.selections.gromacs.SelectionWriter('index_density_sterol_lipid_lower.ndx', mode='w') as ndx:
            ndx.write(self.sterol_head_down, name='lower_sterol_lipid')

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_sterol_lipid_upper, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_sterol_lipid_upper_raw, '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', index_density_sterol_lipid_lower, '-b', str(start_time), '-e', str(end_time), \
            '-o', density_sterol_lipid_lower_raw, '-dens', 'number', '-sl', '100']        
        out, err = exec_gromacs(get_density)
        
        with open("density_sterol.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("density_sterol_lipid_upper_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_sterol_lipid_upper.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
        with open("density_sterol_lipid_lower_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("density_sterol_lipid_lower.dat" , 'w') as fout:
            
            fout.write('Z\tDensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
