'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from . import order
from .order import Order
from .order import tilt_sterol
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
import MDAnalysis as mda
import numpy as np


LOGGER = log.LOGGER

class Tilt_sterol(Order):
    '''
        This class handles the calculation of the tilt of sterol molecules
    '''
    LOGGER = LOGGER

    def __init__(self, inputfilename="inputfile"):
        super().__init__(inputfilename)
    
        self.lipid_type_items = ' '.join(self.molecules)
        self.u = mda.Universe(self.gropath,self.trjpath)  
        self.sterol_lipid = ''.join(self.molecules[1])
                            
    def tilt_calculation(self,start_time,end_time):
        
        n_frame = 0
        with open("tilt_sterol.dat" , 'w') as fout:
                
                fout.write('Time\ttilt_angle\n')
                
                for i_ts,ts in enumerate(self.u.trajectory[start_time:end_time:]):
                    time = self.u.trajectory.time
                    n_frame += 1
                    sterol = self.u.select_atoms('resname {}'.format(self.sterol_lipid))
                    calc_s = tilt_sterol
                    resids_list = list(set(sterol.resids))                    
                    tilt_i = 0
                    for res in resids_list:
                        tilt_value = calc_s(self.u, res)
                        #print(tilt_value)
                        tilt_i += tilt_value
                    tilt_mean = tilt_i / len(resids_list)
                    fout.write('{}\t{}\n'.format(time,tilt_mean))

        