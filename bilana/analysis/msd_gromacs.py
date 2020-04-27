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
            
            self.u = mda.Universe(self.gropath,self.trjpath)
    
    # def MSD_mdanalysis(self,start_frame,end_frame):
    #     '''This function calulate the MSD through MDAnalaysis'''

    #     u = mda.Universe(self.gropath,self.trjpath)
        
    #     selection = ('resname {} and name P'.format(self.lipid_types_mainlipid))
    #     print(self.lipid_types_mainlipid)
        
    #     if self.times[2] == '1000':
    #         start_frame = 100
    #         frame_intervals = 1
    #         end_frame = 300
    #     else:
    #         start_frame = 1000
    #         frame_intervals = 10 
    #         end_frame = 3000       

    #     MSD_analysis = MSD(u, selection, start_frame, end_frame, 20)
    #     MSD_analysis.run()

    #     with open("msd_mdanalysis.xvg" , 'w') as fout:
    #         time = 0
    #         fout.write('Time\tMSD\n')
    #         for msd in MSD_analysis.timeseries:
    #             fout.write("{time} {msd}\n".format(time=time, msd=msd))
    #             time += 1
        
    def MSD_gromacs(self, selection, outputfile_name, start_time, end_time):

        '''This function calulate the MSD through Gromacs'''
        
        sel = self.u.select_atoms('{}'.format(selection))

        with mda.selections.gromacs.SelectionWriter('index_msd.ndx', mode='w') as ndx:
            ndx.write(sel, name='{}'.format(selection))

        get_msd = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', 'index_msd.ndx', \
            '-o', "msd_" + outputfile_name + "_raw", '-lateral', 'z', '-b', str(start_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1']
            
        #print(get_msd)  
        out, err = exec_gromacs(get_msd)

        with open("gmx_msd_" + outputfile_name + ".log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("msd_" + outputfile_name + "_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("msd_" + outputfile_name + ".dat" , 'w') as fout:
            
            fout.write('Time\tMSD\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
    
    def MSD_gromacs_slices(self, selection, outputfile_name, start, end, interval):

        '''This function calulate the MSD through Gromacs'''
        
        sel = self.u.select_atoms('{}'.format(selection))

        with mda.selections.gromacs.SelectionWriter('index_msd.ndx', mode='w') as ndx:
            ndx.write(sel, name='{}'.format(selection))

        for i_n, i in enumerate(np.arange(start,end,interval)):

            start_time, end_time = i, i + interval

            get_msd = [GMXNAME, 'msd', '-f', self.trjpath, '-s', self.tprpath, '-n', 'index_msd.ndx', \
                '-o', "msd_" + outputfile_name + "_raw_{}".format(str(i_n)), '-lateral', 'z', '-b', str(start_time), '-e', str(end_time), '-rmcomm', '-beginfit', '-1', '-endfit', '-1']
                
            #print(get_msd)  
            out, err = exec_gromacs(get_msd)

            with open("gmx_msd_" + outputfile_name + ".log","a") as logfile:
                logfile.write(err)
                logfile.write(out)

            with open("msd_" + outputfile_name + "_raw_{}.xvg".format(str(i_n)), 'r') as f:
                ls = f.readlines()

            with open("msd_" + outputfile_name + "_{}.dat".format(str(i_n)) , 'w') as fout:
                
                fout.write('Time\tMSD\n')
                for l in ls:
                    lc = l.strip()
                    if lc:
                        if lc[0] != '#' and lc[0] != '@':
                            lf = lc.split()
                            fout.write('{}\t{}\n'.format(lf[0],lf[1]))