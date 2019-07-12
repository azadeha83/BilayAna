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
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD
from numpy import linalg

LOGGER = log.LOGGER

class MSDanalysisDirect(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)
    
            self.lipid_type_items = ' '.join(self.molecules)
            self.lipid_types_first = ''.join(self.molecules[0])
    
    def MSD(self,start_frame,end_frame,ref_time):
        '''This function calulate the MSD through MDAnalaysis'''
                
        u = mda.Universe(self.gropath,self.trjpath)

        selection = u.select_atoms('resname {} and name P'.format(self.lipid_types_first))

        for ts in u.trajectory[ref_time]:       
                
            selection_ref_xyz = selection.positions
        print(selection_ref_xyz)

        t0 = ref_time*100
        with open("msd_direct_mdanalysis.xvg" , 'w') as fout:
            
            fout.write('Time\tMSD\n')

            for i_ts,ts in enumerate(u.trajectory[start_frame:end_frame:]):
                time = u.trajectory.time
                #print(len(self.u.trajectory))
                #n_frame += 1
                selection_xyz = selection.positions
                #print(selection_xyz)
                displacement = selection_xyz - selection_ref_xyz
                #print(displacement**2)
                mean_squared_displacement = np.sum(displacement**2) / (selection.n_atoms)

                fout.write('{}\t{}\n'.format(u.trajectory.time - t0,mean_squared_displacement))
            
            
        


