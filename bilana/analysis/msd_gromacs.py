'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from . import neighbors
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
from . import neighbors
from .neighbors import Neighbors
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

LOGGER = log.LOGGER

class MSDanalysis(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)
    
            self.lipid_type_items = ' '.join(self.molecules)
            self.lipid_types_mainlipid = ''.join(self.molecules[0])
            #self.lipid_types_sterol = ''.join(self.molecules[1])
            self.u = mda.Universe(self.gropath,self.trjpath)

    
    def MSD_mdanalysis(self,start_frame,end_frame):
        '''This function calulate the MSD through MDAnalaysis'''

        u = mda.Universe(self.gropath,self.trjpath)
        
        selection = ('resname {} and name P'.format(self.lipid_types_mainlipid))
        print(self.lipid_types_mainlipid)
        
        if self.times[2] == '1000':
            start_frame = 100
            frame_intervals = 1
            end_frame = 300
        else:
            start_frame = 1000
            frame_intervals = 10 
            end_frame = 3000       

        MSD_analysis = MSD(u, selection, start_frame, end_frame, 20)
        MSD_analysis.run()

        with open("msd_mdanalysis.xvg" , 'w') as fout:
            time = 0
            fout.write('Time\tMSD\n')
            for msd in MSD_analysis.timeseries:
                fout.write("{time} {msd}\n".format(time=time, msd=msd))
                time += 1

    def diffusion_gromacs(self, selection, start_time, end_time):

        '''This function calulate the MSD and diffusion coefficient through Gromacs'''

        cwd = os.getcwd()
        sel = self.u.select_atoms('{}'.format(selection))
       
        msd_out = ''.join(cwd + '/' + 'msd' + '_' + '{}'.format('_'.join(map(str, selection.split(' ')))))
        diff_out = ''.join(cwd + '/' + 'diff' + '_' + '{}'.format('_'.join(map(str, selection.split(' ')))))
        
        with mda.selections.gromacs.SelectionWriter('index_msd.ndx', mode='w') as ndx:
            ndx.write(sel, name='{}'.format(selection))

        get_msd = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', 'index_msd.ndx', \
            '-o', msd_out, '-lateral', 'z', '-b', str(start_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1', '-xvg', 'none']
            
        out, err = exec_gromacs(get_msd)

        # get_diff = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', 'index_msd.ndx', \
        #     '-o', msd_out, '-mol', diff_out, '-lateral', 'z', '-b', str(start_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1', '-xvg', 'none']
        # inpstr = '{}'.format(selection) + "\n"
        
        # #get_diff = get_msd + ["-mol", diff_out]
        # out, err = exec_gromacs(get_diff,inpstr)


    def MSD_gromacs_mainlipid(self,start_time,end_time):

        '''This function calulate the MSD through Gromacs'''

        cwd = os.getcwd()
        index_msd = ''.join(cwd + '/' + 'index_msd_mainlipid.ndx')
       
        msd_raw = ''.join(cwd + '/' + 'msd_mainlipid_raw')

        get_selection = [GMXNAME, 'select', '-f', self.gropath, '-s', self.tprpath, '-on', index_msd, \
            '-select', '(resname {} and name P)'.format(self.lipid_types_mainlipid)]

        #print(get_selection)
        
        out, err = exec_gromacs(get_selection)

        get_msd = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', index_msd, \
            '-o', msd_raw, '-lateral', 'z', '-b', str(start_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1']
            
        #print(get_msd)  
        out, err = exec_gromacs(get_msd)

        with open("gmx_msd_mainlipid.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("msd_mainlipid_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("msd_mainlipid.dat" , 'w') as fout:
            
            fout.write('Time\tMSD\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))

    def MSD_gromacs_sterol(self,start_time,end_time):

        '''This function calulate the MSD through Gromacs'''

        cwd = os.getcwd()
        index_msd = ''.join(cwd + '/' + 'index_msd_sterol.ndx')
       
        msd_raw = ''.join(cwd + '/' + 'msd_sterol_raw')
        
        head_atoms = lipidmolecules.head_atoms_of(self.lipid_types_sterol)
        head_atom = head_atoms[0]

        get_selection = [GMXNAME, 'select', '-f', self.gropath, '-s', self.tprpath, '-on', index_msd, \
            '-select', '(resname {} and name {})'.format(self.lipid_types_sterol,head_atom)]

        #print(get_selection)
        
        out, err = exec_gromacs(get_selection)

        get_msd = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', index_msd, \
            '-o', msd_raw, '-lateral', 'z', '-b', str(start_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1']
            
        #print(get_msd)  
        out, err = exec_gromacs(get_msd)

        with open("gmx_msd_sterol.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("msd_sterol_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("msd_sterol.dat" , 'w') as fout:
            
            fout.write('Time\tMSD\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
                        
