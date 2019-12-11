'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from . import order
from .order import Order
from . import neighbors
from .neighbors import Neighbors
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
import MDAnalysis as mda
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftfreq
from matplotlib.backends.backend_pdf import PdfPages
from ..files.io import *

LOGGER = log.LOGGER

class Density(Neighbors):

    def __init__(self,inputfilename="inputfile"):
        super().__init__(inputfilename)
    
        # self.lipid_type_items = ' '.join(self.molecules)
        self.u = mda.Universe(self.gropath,self.trjpath)
        self.neiblist = neighbors.get_neighbor_dict()
        
        # self.main_lipid = ''.join(self.molecules[0])
        # self.sterol_lipid = ''.join(self.molecules[1])
        
        # self.lipid_type_items = ' '.join(self.molecules)
        
        # self.lipid_head = self.u.select_atoms('resname {} and name P'.format(self.main_lipid)) 
        
        # head_atoms = lipidmolecules.head_atoms_of(self.sterol_lipid)
        # head_atom = head_atoms[0]

        # self.sterol_head = self.u.select_atoms('resname {} and name {}'.format(self.sterol_lipid,head_atom)) 
    '''    
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
    '''
    def density(self, selection, outputfile_name, start_time, end_time):
        
        sel = self.u.select_atoms('{}'.format(selection))

        with mda.selections.gromacs.SelectionWriter('index_density.ndx', mode='w') as ndx:
            ndx.write(sel, name='{}'.format(selection))

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-n', 'index_density.ndx', '-b', str(start_time), '-e', str(end_time), \
                    '-o', outputfile_name, '-center', '-dens', 'number', '-sl', '100', '-xvg', 'none']        
        
        out, err = exec_gromacs(get_density)
    
    
    
    def density_2D(self, sterol, sterol_resid, lipid, bin_width, start_frame, end_frame, time_interval):
        ''' This function calculated the density of lipid tails around a specific sterol'''
        
        ref = mda.Universe('{}'.format(self.gropath))         
        print(len(self.u.trajectory))   
        
        headatms = lipidmolecules.head_atoms_of(sterol)[:18]
        print(headatms)      
        
        sterol_sel = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))
        sterol_sel_ref = ref.select_atoms('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))

        print('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))
        print(sterol_sel_ref)

        sterol_ref_center = sterol_sel_ref.center_of_geometry()
        sterol_ref_center -= np.array(sterol_ref_center)
            
        x_min, x_max = sterol_ref_center[0] - 20 , sterol_ref_center[0] + 20
        y_min, y_max = sterol_ref_center[1] - 20 , sterol_ref_center[1] + 20

        print(x_min,x_max)
        
        x_edges = np.arange(x_min,x_max,bin_width)
        y_edges = np.arange(y_min,y_max,bin_width)

        sterol_ref_center = sterol_sel_ref.center_of_geometry()
                       
        n_frame = 0
        H_sterol_edges = 0
        H_sterol_methyl = 0
        H_lipidchains = 0
        volume = bin_width**2
        
        for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:time_interval]):
            time = self.u.trajectory.time
            
            mda.analysis.align.alignto(self.u.atoms, ref.atoms, select='resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))       
            
            sterol_edgeatoms = lipidmolecules.head_atoms_of(sterol)[7] + ' ' + lipidmolecules.head_atoms_of(sterol)[11] 
            sterol_methylatom = lipidmolecules.head_atoms_of(sterol)[18]
            
            C_sterol_edgeatoms = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol, sterol_resid, sterol_edgeatoms)).positions
            C_sterol_edgeatoms -= np.array(sterol_ref_center)

            C_sterol_methylatom = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol, sterol_resid, sterol_methylatom)).positions
            C_sterol_methylatom -= np.array(sterol_ref_center)
            
            hist_edges,x_bins,y_bins = np.histogram2d(C_sterol_edgeatoms[:,0].flatten(),C_sterol_edgeatoms[:,1].flatten(), (x_edges,y_edges))
            H_sterol_edges += hist_edges

            hist_methyl,x_bins,y_bins = np.histogram2d(C_sterol_methylatom[:,0].flatten(),C_sterol_methylatom[:,1].flatten(), (x_edges,y_edges))
            H_sterol_methyl += hist_methyl

            neibs = self.neiblist[sterol_resid][float(time)]
            neib_list = []
            if len(neibs) != 0:
                for res in neibs:
                    resn = self.resid_to_lipid[res]
                    if resn not in lipidmolecules.STEROLS:
                        neib_list.append(res)
            
            flattened_carbons = [y for x in lipidmolecules.tailcarbons_of(lipid) for y in x]
            tail_carbons_xyz = self.u.select_atoms("resname {} and resid {} and name {}".format(lipid, ' '.join(map(str, neib_list)), ' '.join(map(str, flattened_carbons)))).positions
            #print(tail_carbons_xyz)
            tail_carbons_xyz -= np.array(sterol_ref_center)
            hist1,x_bins1,y_bins1 = np.histogram2d(tail_carbons_xyz[:,0].flatten(),tail_carbons_xyz[:,1].flatten(), (x_edges,y_edges))
            H_lipidchains += hist1
            
            n_frame += 1

        H_sterol_edges = np.true_divide(H_sterol_edges,volume)
        H_sterol_methyl = np.true_divide(H_sterol_methyl,volume)
        H_lipidchains = np.true_divide(H_lipidchains,volume)
        H_sterol_edges /= n_frame
        H_sterol_methyl /= n_frame
        H_lipidchains /= n_frame
        
        a = np.array(H_sterol_edges)
        b = np.array(H_sterol_methyl)
        c = np.array(H_lipidchains)
        
        np.save('x_edges_density.dat',x_edges)
        np.save('y_edges_density.dat',y_edges)
        np.save('density_2d_sterol_edges.dat',a)
        np.save('density_2d_sterol_methyl.dat',b)
        np.save('density_2d_lipidchains.dat',c)

    def density_2D_all_sterols(self, sterol, sterol_resid, lipid, bin_width, start_frame, end_frame, time_interval):
        ''' This function calculated the density of lipid tails around a specific sterol'''
        
        ref = mda.Universe('{}'.format(self.gropath))         
        print(len(self.u.trajectory))   
        
        headatms = lipidmolecules.head_atoms_of(sterol)[:18]
        print(headatms)      
        
        sterol_sel = self.u.select_atoms('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))
        sterol_sel_ref = ref.select_atoms('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))

        print('resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))
        print(sterol_sel_ref)

        sterol_ref_center = sterol_sel_ref.center_of_geometry()
        sterol_ref_center -= np.array(sterol_ref_center)
            
        x_min, x_max = sterol_ref_center[0] - 20 , sterol_ref_center[0] + 20
        y_min, y_max = sterol_ref_center[1] - 20 , sterol_ref_center[1] + 20

        print(x_min,x_max)
        
        x_edges = np.arange(x_min,x_max,bin_width)
        y_edges = np.arange(y_min,y_max,bin_width)

        sterol_ref_center = sterol_sel_ref.center_of_geometry()
                       
        n_frame = 0
        H_sterol_edges = 0
        H_sterol_methyl = 0
        H_lipidchains = 0
        volume = bin_width**2
        
        for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:time_interval]):
            time = self.u.trajectory.time

            for res_s in self.MOLRANGE:
                resn_s = self.resid_to_lipid[res_s]
                    
                if resn_s not in lipidmolecules.STEROLS:
                    continue  
                
                mda.analysis.align.alignto(self.u.atoms, ref.atoms, select='resname {} and resid {} and name {}'.format(sterol, sterol_resid, ' '.join(map(str, headatms))))       
                
                sterol_edgeatoms = lipidmolecules.head_atoms_of(sterol)[7] + ' ' + lipidmolecules.head_atoms_of(sterol)[11] 
                sterol_methylatom = lipidmolecules.head_atoms_of(sterol)[18]
                
                C_sterol_edgeatoms = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol, res_s, sterol_edgeatoms)).positions
                C_sterol_edgeatoms -= np.array(sterol_ref_center)

                C_sterol_methylatom = self.u.select_atoms("resname {} and resid {} and name {}".format(sterol, res_s, sterol_methylatom)).positions
                C_sterol_methylatom -= np.array(sterol_ref_center)
                
                hist_edges,x_bins,y_bins = np.histogram2d(C_sterol_edgeatoms[:,0].flatten(),C_sterol_edgeatoms[:,1].flatten(), (x_edges,y_edges))
                H_sterol_edges += hist_edges

                hist_methyl,x_bins,y_bins = np.histogram2d(C_sterol_methylatom[:,0].flatten(),C_sterol_methylatom[:,1].flatten(), (x_edges,y_edges))
                H_sterol_methyl += hist_methyl

                neibs = self.neiblist[res_s][float(time)]
                neib_list = []
                if len(neibs) != 0:
                    for res in neibs:
                        resn = self.resid_to_lipid[res]
                        if resn not in lipidmolecules.STEROLS:
                            neib_list.append(res)
                
                flattened_carbons = [y for x in lipidmolecules.tailcarbons_of(lipid) for y in x]
                tail_carbons_xyz = self.u.select_atoms("resname {} and resid {} and name {}".format(lipid, ' '.join(map(str, neib_list)), ' '.join(map(str, flattened_carbons)))).positions
                #print(tail_carbons_xyz)
                tail_carbons_xyz -= np.array(sterol_ref_center)
                hist1,x_bins1,y_bins1 = np.histogram2d(tail_carbons_xyz[:,0].flatten(),tail_carbons_xyz[:,1].flatten(), (x_edges,y_edges))
                H_lipidchains += hist1
                
                n_frame += 1

        H_sterol_edges = np.true_divide(H_sterol_edges,volume)
        H_sterol_methyl = np.true_divide(H_sterol_methyl,volume)
        H_lipidchains = np.true_divide(H_lipidchains,volume)
        H_sterol_edges /= n_frame
        H_sterol_methyl /= n_frame
        H_lipidchains /= n_frame
        
        a = np.array(H_sterol_edges)
        b = np.array(H_sterol_methyl)
        c = np.array(H_lipidchains)
        
        np.save('x_edges_density.dat',x_edges)
        np.save('y_edges_density.dat',y_edges)
        np.save('density_2d_sterol_edges_all.dat',a)
        np.save('density_2d_sterol_methyl_all.dat',b)
        np.save('density_2d_lipidchains_all.dat',c)
    
    def orient(self, resid, c_sterol_xyz, neib_xyz):
        ''' Calculate the orientation angle of sterol molecule with its neighbors '''
                
        resn = self.resid_to_lipid[resid]

        C8_sterol = self.u.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[8])).positions
        C10_sterol = self.u.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[10])).positions
        C13_sterol = self.u.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[13])).positions
        #C19_sterol = mda_uni.select_atoms("resid {} and name {}".format(resid, lipidmolecules.head_atoms_of(resn)[19])).positions # just for checking
        
        v1 = np.subtract(C8_sterol,C10_sterol)
        v2 = np.subtract(C13_sterol,C10_sterol)
        v3 = np.cross(v1,v2) # vector peripendicular to the plane of sterol molecule
        #print(neib_xyz)
        #print(c_sterol_xyz)
        a = np.subtract(neib_xyz,c_sterol_xyz)
        #print(a)
        a[0][2] = 0

        orient_angle = np.arccos(np.dot(a[0], v3[0])/(np.linalg.norm(a)*np.linalg.norm(v3))) * (180/np.pi)
        
        orient_flag = 0
        
        if orient_angle < 90:
            orient_flag += 1
        
        return orient_flag
    
    def density_packing(self, sterol, lipid, start_frame, end_frame, time_interval):
        ''' This function calculated the density of lipid tails around a specific sterol'''
        
        flattened_carbons = [y for x in lipidmolecules.tailcarbons_of(lipid) for y in x]
        print(lipidmolecules.head_atoms_of(sterol)[1:] + lipidmolecules.tail_atoms_of(sterol)[0])
        carbons_sterol =  list(set(lipidmolecules.head_atoms_of(sterol)[1:] + lipidmolecules.tail_atoms_of(sterol)[0]))
        
        with open('density_packing.dat', "w") as out:
            print("{: <12}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15}".format("Time", "Resid", "Atomname", "N_lipid", "N_sterol", "N_smooth", "N_rough", "Avg_distance"\
                , "Avg_dist_rough", "Avg_dist_smooth") ,file=out)
            
            for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:time_interval]):
                time = self.u.trajectory.time
                print(time)
                
                for res_s in self.MOLRANGE:
                        resn_s = self.resid_to_lipid[res_s]
                        if resn_s not in lipidmolecules.STEROLS:
                            continue            
                        
                        for i_c in carbons_sterol:
                            #print("resname {} and resid {} and name {}".format(resn_s, res_s, ''.join(map(str, i_c))))
                            C_sterol_xyz = self.u.select_atoms("resname {} and resid {} and name {}".format(resn_s, res_s, ''.join(map(str, i_c)))).positions

                            neibs = self.neiblist[res_s][float(time)]
                            
                            smooth = 0
                            rough = 0
                            n_lipid = 0
                            n_sterol = 0
                            sum_dist = 0
                            sum_dist_smooth = 0
                            sum_dist_rough = 0
                            
                            if len(neibs) != 0:
                                for res_n in neibs:
                                    resn_n = self.resid_to_lipid[res_n]
                                    
                                    if resn_n not in lipidmolecules.STEROLS:
                                        
                                        for i in flattened_carbons:
                                            
                                            C_xyz = self.u.select_atoms("resname {} and resid {} and name {}".format(resn_n, res_n, ''.join(map(str, i)))).positions
                                            #print("resname {} and resid {} and name {}".format(resn_n, res_n, ''.join(map(str, i))))
                                            
                                            dist = np.linalg.norm((C_sterol_xyz - C_xyz)[0:3])
                                            
                                            if dist < 7.0:
                                                
                                                sum_dist += dist
                                                n_lipid += 1
                                                
                                                if self.orient(res_s, C_sterol_xyz, C_xyz) == 1:
                                                    rough += 1
                                                    sum_dist_rough += dist
                                                else:
                                                    smooth += 1
                                                    sum_dist_smooth += dist

                                    elif resn_n in lipidmolecules.STEROLS:
                                        for i in carbons_sterol:
                                            
                                            C_s_xyz = self.u.select_atoms("resname {} and resid {} and name {}".format(resn_n, res_n, ''.join(map(str, i)))).positions
                                            #print("resname {} and resid {} and name {}".format(resn_n, res_n, ''.join(map(str, i))))
                                            
                                            dist = np.linalg.norm((C_sterol_xyz - C_s_xyz)[0:3])
                                            
                                            if dist < 7.0:
                                                
                                                sum_dist += dist                                           
                                                n_sterol += 1
                                                
                                                if self.orient(res_s, C_sterol_xyz, C_s_xyz) == 1:
                                                    rough += 1
                                                    sum_dist_rough += dist
                                                else:
                                                    smooth += 1
                                                    sum_dist_smooth += dist

                                if (n_sterol + n_lipid) != 0:
                                    avg_dist = sum_dist / (n_sterol + n_lipid)  
                                else:
                                    avg_dist = np.nan                          
                                
                                if rough !=0:
                                    avg_dist_rough = sum_dist_rough / rough
                                else:
                                     avg_dist_rough = np.nan  

                                if smooth !=0:
                                    avg_dist_smooth = sum_dist_smooth / smooth
                                else:
                                    avg_dist_smooth = np.nan
                            else:
                                avg_dist = np.nan                            
                                avg_dist_rough = np.nan                            
                                avg_dist_smooth = np.nan

                                
                            print("{: <12.2f}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15}{: <15.2f}{: <15.2f}{: <15.2f}".format(time, res_s, ''.join(map(str, i_c)), n_lipid, n_sterol, smooth, rough, avg_dist,\
                                avg_dist_rough, avg_dist_smooth), file=out)

    
    def density_packing_molecules(self, sterol, lipid, start_frame, end_frame, time_interval):
        ''' This function calculated the density of lipid around a specific sterol'''
        calc_s = Order.scc_of_res
        
        with open('density_packing_molecules.dat', "w") as out, open('Number_of_DPPC_neibors.dat', "w") as out1:
            print("{: <12}{: <15}{: <15}{: <15}{: <15}{: <15}".format("Time", "Resid", "Neib_type", "Orient_flag", "Distance", "Scd"\
                ) ,file=out)
            print("{: <12}{: <15}{: <15}{: <15}{: <15}{: <15}".format("Time", "N-smooth", "N_rough", "Scd-smooth", "Scd_rough", "N_sterol"), file=out1)
            
            for i_ts,ts in enumerate(self.u.trajectory[start_frame:end_frame:time_interval]):
                time = self.u.trajectory.time
                print(time)

                smooth = 0
                rough = 0
                n_sterol = 0
                scd_rough = 0
                scd_smooth = 0

                n_lipid_smooth = 0
                n_lipid_rough = 0
                
                for res_s in self.MOLRANGE:
                    resn_s = self.resid_to_lipid[res_s]
                    
                    if resn_s not in lipidmolecules.STEROLS:
                        continue            
                    
                        #print("resname {} and resid {} and name {}".format(resn_s, res_s, ''.join(map(str, i_c))))
                    O_sterol_xyz = self.u.select_atoms("resname {} and resid {} and name O3".format(resn_s, res_s)).positions

                    neibs = self.neiblist[res_s][float(time)]                        
                    
                    orient_flag = 0
                    distance = 0 


                    if len(neibs) != 0:
                        
                        n_sterol += 1
                                               
                        for res_n in neibs:
                            resn_n = self.resid_to_lipid[res_n]

                            if resn_n in lipidmolecules.STEROLS:
                                continue
                                                               
                            C_xyz = self.u.select_atoms("resname {} and resid {} and name P".format(resn_n, res_n,)).positions
                            #print("resname {} and resid {} and name {}".format(resn_n, res_n, ''.join(map(str, i))))
                            
                            dist = np.linalg.norm((O_sterol_xyz - C_xyz)[0:2])
                            scd_value = calc_s(self.u, int(res_n))
                            
                            if self.orient(res_s, O_sterol_xyz, C_xyz) == 1:
                                rough += 1
                                orient_flag = 1
                                scd_rough += scd_value
                                n_lipid_rough += 1
                                distance = dist
                            else:
                                smooth += 1
                                orient_flag = 0
                                scd_smooth += scd_value
                                n_lipid_smooth += 1
                                distance = dist 
                            
                            
                            print("{: <12.2f}{: <15}{: <15}{: <15}{: <15.2f}{: <15.2f}".format(time, res_s, resn_n, orient_flag, distance, scd_value,\
                                ), file=out)
                    
                print("{: <12.2f}{: <15.2f}{: <15.2f}{: <15.2f}{: <15.2f}{: <15}".format(time, smooth/n_sterol, rough/n_sterol, scd_smooth/n_lipid_smooth, scd_rough/n_lipid_rough, n_sterol), file=out1)                
                            