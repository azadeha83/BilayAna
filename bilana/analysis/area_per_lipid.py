'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import os
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
import MDAnalysis as mda
import pandas as pd
import numpy as np
import csv

LOGGER = log.LOGGER

class Areaperlipid_compressibility(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)

            self.lipid_type_items = ' '.join(self.molecules)
            self.lipid_types_first = ''.join(self.molecules[0])

    def area_per_lipid(self,end_time):
        ''' This module calculates the area per lipid considering all the lipid molecules including sterol molecules as well as compresseibility'''

        cwd = os.getcwd()
        edrout = ''.join(cwd + '/' + 'box_xy')

        box_size_arglist = [GMXNAME, 'energy', '-f', self.edrpath, '-o', edrout, '-e', str(end_time)]

        inp_str = b'Box-X\nBox-Y\n0\n'
        out, err = exec_gromacs(box_size_arglist, inp_str)

        with open("gmx_box_size.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        number_of_lipids = self.number_of_lipids/2
        print(number_of_lipids)

        with open("box_xy.xvg", 'r') as f:
            ls = f.readlines()

        with open("area_per_lipid.dat" , 'w') as fout:

            fout.write('Time\tbox_x\tbox_y\tarea_per_lipid\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        area_per_lipid = (float(lf[1])*float(lf[2]))/number_of_lipids
                        fout.write('{}\t{}\t{}\t{}\n'.format(lf[0],lf[1],lf[2],area_per_lipid))

    def area_per_lipid_mainlipid(self,end_time):

        ''' This module calculates the area per lipid considering only the main lipid molecule'''

        u = mda.Universe(self.tprpath,self.trjpath)

        main_lipid_selection = u.select_atoms('resname {} and name P'.format(self.lipid_types_first))
        number_of_main_lipid_selection = len(main_lipid_selection)
        print(number_of_main_lipid_selection)
        LOGGER.debug("number of main lipids: %d", number_of_main_lipid_selection)

        cwd = os.getcwd()
        edrout = ''.join(cwd + '/' + 'box_xy')

        box_size_arglist = [GMXNAME, 'energy', '-f', self.edrpath, '-o', edrout, '-e', str(end_time)]
        
        inp_str = b'Box-X\nBox-Y\n0\n'
        out, err = exec_gromacs(box_size_arglist, inp_str)

        with open("gmx_box_size.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        number_of_lipids = number_of_main_lipid_selection/2
        print(number_of_lipids)

        with open("box_xy.xvg", 'r') as f:
            ls = f.readlines()

        with open("area_per_lipid_mainlipid.dat" , 'w') as fout:

            fout.write('Time\tbox_x\tbox_y\tarea_per_lipid\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        area_per_lipid = (float(lf[1])*float(lf[2]))/number_of_lipids
                        fout.write('{}\t{}\t{}\t{}\n'.format(lf[0],lf[1],lf[2],area_per_lipid))

    def compressibility(self,start_time,number_of_blocks):

        ''' This module calculates the using the are_per_lipid data calculated in the first function of this class'''

        df = pd.read_table('area_per_lipid.dat', header=0, delim_whitespace=True)
        df1 = df.loc[df['Time'] > start_time]
        #df1 = df1.loc[df1['Time'] % 1000 == 0]
        print(len(df1.index))
        n_blocks = number_of_blocks
        area_data = list(df1.loc[:,'area_per_lipid'])
        l = int(len(area_data)/n_blocks)
        block_lists = []
        
        j = 0
        for i in range(n_blocks):
            blocks = area_data[i+j:i+j+l-1]
            block_lists.append(blocks)
            j += l-1

        area_blocks = np.array(block_lists)

        area_mean = np.mean(area_blocks,axis=1)
        area_variance = np.var(area_blocks,axis=1)
        
        T = float(self.temperature)
        Boltzman = 1.380649 * 1e-23 
        comp = Boltzman*T*(area_mean/area_variance) * 1e+18

        with open('area_compressibility.csv','w') as f:
                    writer = csv.writer(f,delimiter='\t')
                    writer.writerow(["mean_area", "std_area", "compressibility", "std_compressibility", "se_compressibility"])
                    writer.writerow([np.mean(area_mean),np.std(area_mean),np.mean(comp),np.std(comp),np.std(comp)/np.sqrt(n_blocks)])

        # self.area_mean = np.mean(list(df1.iloc[:,1]))
        # self.area_variance = np.var(list(df1.iloc[:,1]))
        # T = float(self.temperature)
        # Boltzman = (1.380649 / 100)*1e-3 # in order to change the units of KT to pN.nm

        # for i,row in df.iterrows():
        #     area = df.at[i,'area_per_lipid']
        #     print(area)
        #     comp = Boltzman*T*(area/self.area_variance)
        #     print(comp)

        #     df.at[i,'compressibility'] = comp
        # df.to_csv('areaperlipid_compress.dat', sep='\t', index=None)
