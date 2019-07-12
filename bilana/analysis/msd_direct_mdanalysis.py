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

        u.trajectory[ref_time]
        selection_ref_xyz = selection.positions
        print(selection_ref_xyz)
        # for ts in u.trajectory[ref_time]:       
                
        #     selection_ref_xyz = selection.positions
        # print(selection_ref_xyz)

        t0 = ref_time*100
        with open("msd_direct_mdanalysis.xvg" , 'w') as fout:
            
            fout.write('Time\tMSD\n')

            for i_ts,ts in enumerate(u.trajectory[ref_time:end_frame:]):
                time = u.trajectory.time
                #print(len(self.u.trajectory))
                #n_frame += 1
                selection_xyz = selection.positions
                print(selection_xyz)
                displacement = selection_xyz - selection_ref_xyz
                print(displacement)

                displacement_norm = np.linalg.norm(displacement, axis=1)
                print(displacement_norm)
                mean_squared_displacement = np.mean(1e-2*(displacement_norm**2))
                print(mean_squared_displacement)

                fout.write('{}\t{}\n'.format(u.trajectory.time - t0,mean_squared_displacement))

    def MSD_smooth(self,start_frame,end_frame,n_ref):
        '''This function calulate the MSD through MDAnalaysis'''
                
        u = mda.Universe(self.gropath,self.trjpath)

        selection = u.select_atoms('resname {} and name P'.format(self.lipid_types_first))
        
        data = np.zeros((len(u.trajectory[start_frame:end_frame:]),n_ref)) 
        time = np.zeros(len(u.trajectory[start_frame:end_frame:]))

        for i in range(n_ref): 
            
            #t0 = (start_frame+i)*100
            t0 = start_frame+i

            u.trajectory[start_frame+i] 

            selection_ref_xyz = selection.positions

            print(selection_ref_xyz)
            n_frame = 0

            for i_ts,ts in enumerate(u.trajectory[start_frame+i:end_frame:]):
                #time = u.trajectory.time
                print(selection.n_atoms)
                #print(len(self.u.trajectory))
                selection_xyz = selection.positions
                # print(selection_xyz)
                # print(selection_ref_xyz)
                displacement = selection_xyz - selection_ref_xyz
                #print(displacement)
                displacement_norm = np.linalg.norm(displacement, axis=1)
                #print(displacement_norm)
                mean_squared_displacement = np.mean(1e-2*displacement_norm)
                #print(mean_squared_displacement)
                #print(int(((time+i) - t0)/10000))
                #data[int(((time+i) - t0)/10000)][i] = mean_squared_displacement
                time[start_frame+i+n_frame-t0] = n_frame*100
                data[start_frame+i+n_frame-t0][i] = mean_squared_displacement
                #print(data)
                n_frame += 1

        
        final_msd = np.true_divide(data.sum(1),(data!=0).sum(1))
        total_data = np.vstack((time,final_msd))
        print(total_data)

        with open("msd_direct_smooth.dat" , 'w') as fout:
            
            fout.write('Time\tMSD\n')
            for i,j in enumerate(total_data[0,:]):
                fout.write('{}\t{}\n'.format(total_data[0,i],total_data[1,i]))
            
        
    def MSD_another(self,start_frame,end_frame,ref_time):
        '''This function calulate the MSD through MDAnalaysis'''
                
        u = mda.Universe(self.gropath,self.trjpath)

        selection = u.select_atoms('resname {} and name P'.format(self.lipid_types_first))

        u.trajectory[ref_time]
        selection_ref_xyz = selection.positions
        print(selection_ref_xyz)
        # for ts in u.trajectory[ref_time]:       
                
        #     selection_ref_xyz = selection.positions
        # print(selection_ref_xyz)

        t0 = ref_time*100
        with open("msd_direct_mdanalysis.xvg" , 'w') as fout:
            
            fout.write('Time\tMSD\n')

            for i_ts,ts in enumerate(u.trajectory[ref_time:end_frame:]):
                time = u.trajectory.time
                #print(len(self.u.trajectory))
                #n_frame += 1
                selection_xyz = selection.positions
                print(selection_xyz)
                displacement = selection_xyz - selection_ref_xyz
                print(displacement)

                displacement_norm = np.linalg.norm(displacement, axis=1)
                print(displacement_norm**2)
                mean_squared_displacement = np.mean(1e-2*displacement_norm**2)
                print(mean_squared_displacement)

                fout.write('{}\t{}\n'.format(u.trajectory.time - t0,mean_squared_displacement))


