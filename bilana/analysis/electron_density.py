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

LOGGER = log.LOGGER

class ElectronDensity(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)
    
            self.lipid_type_items = ' '.join(self.molecules)
    
    def electron_density_mdanalysis(self):
        '''This function calulate the charge density through MDAnalaysis'''

        u = mda.Universe(self.tprpath,self.trjpath)

        lipid_types_first = ''.join(self.molecules[0])
        print(lipid_types_first)
        self.lipid_type_items = ' '.join(self.molecules)

        print(lipid_types_first)

        selection = u.select_atoms('resname {}'.format(self.lipid_type_items))
        
        ref_selection = u.select_atoms('resname {} and name P'.format(lipid_types_first))
        
        print(self.times[2])

        interval_of_frames = str(self.times[2])
        print(interval_of_frames)
        
        if self.times[2] == '1000':
            start_frame = 100
            frame_intervals = 1
            end_frame = 300
        else:
            start_frame = 1000
            frame_intervals = 10 
            end_frame = 3000       

        Idens = LinearDensity(selection, grouping='atoms', binsize=0.25, start=start_frame, step=frame_intervals, stop=end_frame)
        Idens.run()
        Idens.save(description='densprof', form='txt')

        # subtracting the center of geometry of the lipid from the first column of the data
        # calculating center of geometry of lipids
        
        u = self.universe
        cog = []
        for i_ts,ts in enumerate(u.trajectory[start_frame::frame_intervals]):
            
            ref_selection_xyz = ref_selection.center_of_geometry()
            cog.append(ref_selection_xyz[2])

        centerofgeometry_z = np.mean(cog)  
        ####

        dat = np.loadtxt('{}_{}.densprof_atoms.ldens'.format(self.system,self.temperature, skiprows=2))
        dat[:,0] = dat[:,0] - centerofgeometry_z
        np.savetxt('{}_{}.densprof_atoms_centered.ldens'.format(self.system,self.temperature), dat, delimiter=' ', fmt='%.5f')

        # Plotting script
        '''
        df = pd.read_table('{}_{}.densprof_atoms_centered.ldens'.format(self.mdfilepath, self.system), skiprows= [0,1], header=None, delim_whitespace=True)
        
        f = plt.figure()
        ax = f.add_subplot(111)
        
        y = df.loc[:,11].rolling(9, center=True).mean()
        
        ax.plot(df.loc[:,0], (-1)*y, '-r')
        ax.plot(df.loc[:,0], -1* df.loc[:,11], '-b')
        
        plt.xlabel('X [nm]', fontsize=16)
        plt.ylabel('Electron Density [e/', fontsize=16)
        plt.show()
        '''
    def electron_density_gromacs(self,start_time,end_time):

        '''This function calulate the charge density through Gromacs'''

        cwd = os.getcwd()
        index_electrondensity = ''.join(cwd + '/' + 'index_electrondensity.ndx')
        electrons = ''.join(cwd + '/' + 'electrons_partial.dat')
        electron_density_raw = ''.join(cwd + '/' + 'electron_density_raw')

        # get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_electrondensity, \
        #     '-select', '(resname {} TIP3)'.format(self.lipid_type_items) ]

        # print(get_selection)

        #out, err = exec_gromacs(get_selection)

        get_density = [GMXNAME, 'density', '-f', self.trjpath, '-s', self.tprpath, '-b', str(start_time), '-e', str(end_time), \
            '-o', electron_density_raw, '-dens', 'electron', '-ei', electrons, '-center', '-relative']
            
        print(get_density)  
        out, err = exec_gromacs(get_density)

        with open("gmx_density.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("electron_density_raw.xvg", 'r') as f:
            ls = f.readlines()

        with open("electron_density.xvg" , 'w') as fout:
            
            fout.write('Zbins\tdensity\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))

    
'''
def electron_density_gromacs(systeminfo,start_time,end_time):
    This function calulate the charge density through Gromacs
    cwd = os.getcwd()
    lipid_types_first = ''.join(systeminfo.molecules[0])
    print(lipid_types_first)
    lipid_type_items = ' '.join(systeminfo.molecules)
    print(lipid_type_items)

    index_electrondensity = ''.join(cwd + '/' + 'index_electrondensity.ndx')
    electrons = ''.join(cwd + '/' + 'electrons.dat')
    electron_density_raw = ''.join(cwd + '/' + 'electron_density_raw.xvg')
    get_selection = [GMXNAME, 'select', '-f', systeminfo.trjpath, '-s', systeminfo.tprpath, '-on', index_electrondensity, \
        '-select', '"(resname {} TIP3)"'.format(lipid_type_items) ]
    out, err = exec_gromacs(get_selection)
    
    get_density = [GMXNAME, 'density', '-f', systeminfo.trjpath, '-s', systeminfo.tprpath, '-n', index_electrondensity, '-b', start_time, 'e', end_time, \
        '-o', electron_density_raw, '-dens', 'electron', '-ei', electrons, 'center', 'relative', '"(resname {} TIP3)"'.format(lipid_type_items) ]
        
    out, err = exec_gromacs(get_density)
    with open("gmx_density.log","a") as logfile:
        logfile.write(err)
        logfile.write(out)
    number_of_lipids = systeminfo.number_of_lipids/2
    print(number_of_lipids)
    with open("electron_density_raw.xvg", 'r') as f:
        ls = f.readlines()
    with open("electron_density.xvg" , 'w') as fout:
        
        fout.write('Zbins\tdensity\n')
        for l in ls:
            lc = l.strip()
            if lc:
                if lc[0] != '#' and lc[0] != '@':
                    lf = lc.split()
                    fout.write('{}\t{}\n'.format(lf[0],lf[1]))
'''
