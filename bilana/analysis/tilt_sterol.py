'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from . import order
from .order import Order
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
        print(self.sterol_lipid)
    '''                        
    def tilt_sterol(self,start_time,end_time):
        
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
    '''
    def tilt_data(self,start_frame,end_frame):
        n_frame = 0
        tilt_total = []
        with open("tilt_sterol_direct.dat" , 'w') as fout:
                
            fout.write('Time\ttilt_angle\n')
            for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:]):
                time = self.u.trajectory.time
                n_frame += 1
                sterol = self.u.select_atoms('resname {}'.format(self.sterol_lipid))
                print(sterol)
                #print(sterol)
                resids_list = list(set(sterol.resids)) 
                print(resids_list)                   
                tilt_i = 0
                for res in resids_list:
                    #print(res)
                    tilt_value = self.tilt_calculation(self.u, res)
                    #print(tilt_value)
                    tilt_i += tilt_value
                    tilt_total.append(tilt_value)
                    #print(tilt_i)
                tilt_mean = tilt_i / len(resids_list)
                #print(tilt_mean)
                fout.write('{}\t{}\n'.format(time,tilt_mean)) 
        np.savetxt('tilt_total.dat',np.array(tilt_total), delimiter=' ')  

    def tilt_calculation(self, mda_uni, resid):
        ''' Calculate the tilt angle of sterol molecule   '''
        
        resinfo = mda_uni.atoms.select_atoms("resid {}".format(resid))
        #print(resinfo)
        resname = list(set(resinfo.resnames))[0]
        #print(resname)
        tailatms = lipidmolecules.scd_tail_atoms_of(resname)
        atm1, atm2 = tailatms[0][0], tailatms[0][1]
        #print(atm1,atm2)
        coords1 = mda_uni.select_atoms("resid {}  and name {}".format(resid,atm1)).positions
        coords2 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm2)).positions
        #print(coords1)
        diffvector = np.subtract(coords2[0],coords1[0])
        tiltangle = np.arccos(np.dot(diffvector, [0,0,1])/(np.linalg.norm(diffvector))) * (180/np.pi)
        if abs(tiltangle) > 90:
            tiltangle = 180 - tiltangle
        #print(tiltangle)
        return tiltangle
