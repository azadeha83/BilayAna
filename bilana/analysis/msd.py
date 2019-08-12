'''
    This module should contain a class that automates gromacs tool gmx msd
'''
import os
import MDAnalysis as mda
import numpy as np

from MDAnalysis.analysis.lineardensity import LinearDensity
from MDAnalysis.analysis.waterdynamics import MeanSquareDisplacement as MSD

from . import neighbors
from .neighbors import Neighbors
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules

class MSDanalysis(SysInfo):

    def __init__(self,inputfilename="inputfile"):
            super().__init__(inputfilename)

            self.lipid_type_items = ' '.join(self.molecules)
            self.lipid_types_mainlipid = ''.join(self.molecules[0])
            self.lipid_types_sterol = ''.join(self.molecules[1])


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

    def MSD_gromacs_mainlipid(self,start_time,end_time):

        '''This function calulate the MSD through Gromacs'''

        cwd = os.getcwd()
        index_msd = ''.join(cwd + '/' + 'index_msd_mainlipid.ndx')

        msd_raw = ''.join(cwd + '/' + 'msd_mainlipid_raw')

        get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_msd, \
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

        get_selection = [GMXNAME, 'select', '-f', self.trjpath, '-s', self.tprpath, '-on', index_msd, \
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
