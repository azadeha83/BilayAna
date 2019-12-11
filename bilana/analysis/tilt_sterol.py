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
    def tilt_data(self,start_frame,end_frame,step):
        n_frame = 0
        ptilt_total = []
        ttilt_total = []
        with open("tilt_sterol_direct.dat" , 'w') as fout, open('tilt_total.dat', 'w') as fout1:
                
            fout.write('Time\tplane_tilt_angle\ttail_tilt_angle\n')
            fout1.write('Time\tresid\tplane_tilt_angle\ttail_tilt_angle\n')
            
            for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
                time = self.u.trajectory.time
                n_frame += 1
                sterol = self.u.select_atoms('resname {}'.format(self.sterol_lipid))
                resids_list = list(set(sterol.resids)) 
                ptilt_i = 0
                ttilt_i = 0
                
                for res in resids_list:
                    ptilt_value, ttilt_value = self.tilt_calculation(self.u, res)
                    ptilt_i += ptilt_value
                    ttilt_i += ttilt_value
                    ptilt_total.append(ptilt_value)
                    ttilt_total.append(ttilt_value)
                    fout1.write('{: ^8}\t{: ^8}\t{: ^8}\t{: ^8}\n'.format(time, res, ptilt_value, ttilt_value)) 
                
                ptilt_mean = ptilt_i / len(resids_list)
                ttilt_mean = ttilt_i / len(resids_list)
                fout.write('{: ^8}\t{: ^8}\t{: ^8}\n'.format(time,ptilt_mean,ttilt_mean)) 
        
        # np.savetxt('ptilt_total.dat',np.array(ptilt_total), delimiter=' ')  
        # np.savetxt('ttilt_total.dat',np.array(ttilt_total), delimiter=' ')  

    def tilt_calculation(self, mda_uni, resid):
        ''' Calculate the tilt angle of sterol molecule   '''
        
        resinfo = mda_uni.atoms.select_atoms("resid {}".format(resid))

        resname = list(set(resinfo.resnames))[0]

        tailatms1, tailatms2 = lipidmolecules.scd_tail_atoms_of(resname), lipidmolecules.tail_atoms_of(resname)
        atm1, atm2 = tailatms1[0][0], tailatms1[0][1]
        atm3, atm4 = tailatms2[0][5], tailatms2[0][9]


        coords1 = mda_uni.select_atoms("resid {}  and name {}".format(resid,atm1)).positions
        coords2 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm2)).positions

        coords3 = mda_uni.select_atoms("resid {}  and name {}".format(resid,atm3)).positions
        coords4 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm4)).positions
        
        diffvector1 = np.subtract(coords2[0],coords1[0])
        diffvector2 = np.subtract(coords4[0],coords3[0])
        
        ptiltangle = np.arccos(np.dot(diffvector1, [0,0,1])/(np.linalg.norm(diffvector1))) * (180/np.pi)
        ttiltangle = np.arccos(np.dot(diffvector2, [0,0,1])/(np.linalg.norm(diffvector2))) * (180/np.pi)
        
        if abs(ptiltangle) > 90:
            ptiltangle = 180 - ptiltangle
        if abs(ttiltangle) > 90:
            ttiltangle = 180 - ttiltangle

        return ptiltangle, ttiltangle
    
    def metylangle_data(self,start_frame,end_frame,step):
        n_frame = 0
        metylangle_total = []
        with open("metylangle_sterol.dat" , 'w') as fout:
                
            fout.write('Time\tmetylangle\n')
            
            for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:step]):
                time = self.u.trajectory.time
                n_frame += 1
                sterol = self.u.select_atoms('resname {}'.format(self.sterol_lipid))
                
                resids_list = list(set(sterol.resids)) 
                metylangle_i = 0
                
                for res in resids_list:
                    metylangle_value = self.metylangle_calc(self.u, res)
                    metylangle_i += metylangle_value
                    metylangle_total.append(metylangle_value)
                
                metylangle_mean = metylangle_i / len(resids_list)
                
                fout.write('{}\t{}\n'.format(time,metylangle_mean)) 
        
        np.savetxt('metylangle_total.dat',np.array(metylangle_total), delimiter=' ')  
    
    def metylangle_calc(self, mda_uni, resid):
        ''' Calculate the angle between the two metyl groups of the sterol molecule projected on xy plane  '''
        
        resinfo = mda_uni.atoms.select_atoms("resid {}".format(resid))
        #print(resinfo)
        resname = list(set(resinfo.resnames))[0]
        #print(resname)
        planeatms = lipidmolecules.head_atoms_of(resname)
        
        atm1, atm2, atm3, atm4 = planeatms[7], planeatms[10], planeatms[8], planeatms[9]
        
        coords1 = mda_uni.select_atoms("resid {} and name {}".format(resid,atm1)).positions
        coords2 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm2)).positions
        coords3 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm3)).positions
        coords4 = mda_uni.select_atoms("resid {} and name {} ".format(resid,atm4)).positions
        
        diffvector1 = np.multiply(np.subtract(coords2[0],coords1[0]),[1,1,0])
        diffvector2 = np.multiply(np.subtract(coords4[0],coords3[0]),[1,1,0])
        
        metylangle = np.arccos(np.dot(diffvector1, diffvector2)/(np.linalg.norm(diffvector1)*np.linalg.norm(diffvector2))) * (180/np.pi)
        
        return metylangle
