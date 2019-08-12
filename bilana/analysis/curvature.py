'''

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
from scipy.fftpack import fft, fftfreq
from matplotlib.backends.backend_pdf import PdfPages

LOGGER = log.LOGGER

class Curvature(SysInfo):

    def __init__(self,bin_width,inputfilename="inputfile"):
        super().__init__(inputfilename)

        self.lipid_type_items = ' '.join(self.molecules)
        self.bin_width = bin_width
        self.u = mda.Universe(self.gropath,self.trjpath)
        #print(len(self.u.trajectory))

        w = mda.Universe(self.gropath)

        main_lipid = ''.join(self.molecules[0])
        sterol_lipid = ''.join(self.molecules[1])

        self.lipid_type_items = ' '.join(self.molecules)


        main_lipid_selection = w.select_atoms('resname {} and name P'.format(main_lipid))

        main_lipid_selection_center = main_lipid_selection.center_of_geometry()

        main_lipid_xyz = main_lipid_selection.positions


        self.x_min = min(main_lipid_xyz[:,0])
        self.x_max = max(main_lipid_xyz[:,0]) - 15

        self.lipid_head_up = self.u.select_atoms('resname {} and name P O11 O12 O13 O14 and (prop z > {})'.format(main_lipid,main_lipid_selection_center[2]))
        self.lipid_head_down = self.u.select_atoms('resname {} and name P O11 O12 O13 O14 and (prop z < {})'.format(main_lipid,main_lipid_selection_center[2]))

        self.sterol_selection = self.u.select_atoms('resname {} and name O3'.format(sterol_lipid))

        self.x_edges = np.arange(self.x_min,self.x_max,bin_width)

        self.x_bins = 0.5*(self.x_edges[:-1] + self.x_edges[1:])
        self.x = self.x_edges[:]

    def zero_to_nan(self,values):
        return [float('nan') if x==0 else x for x in values]

    n_frame = 0

    def area_yz(self):

        cwd = os.getcwd()
        edrout = ''.join(cwd + '/' + 'box_yz')

        box_size_arglist = [GMXNAME, 'energy', '-f', self.edrpath, '-o', edrout]

        inp_str = b'Box-Y\nBox-Z\n0\n'
        out, err = exec_gromacs(box_size_arglist, inp_str)

        with open("gmx_box_size_yz.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open("box_yz.xvg", 'r') as f:
            ls = f.readlines()

        with open("area_yz.dat" , 'w') as fout:

            fout.write('Time\tvolume\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        area = (float(lf[1])*float(lf[2]))
                        fout.write('{}\t{}\n'.format(lf[0],area))

        df = pd.read_table('area_yz.dat', header=0, delim_whitespace=True)
        df1 = df.loc[df['Time'] > 100000]
        self.area_mean = np.mean(list(df1.iloc[:,1]))
        return self.area_mean

    def fourier_transform(self,start_time,end_time):

        density_profiles_sterol = []
        z_profiles_up = []
        z_profiles_down = []

        n_frame = 0

        for i_ts,ts in enumerate(self.u.trajectory[start_time:end_time:]):
            time = self.u.trajectory.time
            #print(len(self.u.trajectory))
            n_frame += 1
            head_up_xyz = self.lipid_head_up.positions
            head_down_xyz = self.lipid_head_down.positions

            sterol_xyz = self.sterol_selection.positions

            head_up_x = head_up_xyz[:,0]
            head_down_x = head_down_xyz[:,0]
            head_up_z = head_up_xyz[:,2]
            head_down_z = head_down_xyz[:,2]
            head_sterol_x = sterol_xyz[:,0]

            ib_up = np.rint(np.divide((head_up_x), self.bin_width))
            ib_down = np.rint(np.divide((head_down_x), self.bin_width))
            ib_sterol = np.rint(np.divide((head_sterol_x), self.bin_width))
            bins = int((self.x_max-self.x_min)/self.bin_width)
            #print(len(range(bins)))


            z_up = np.empty(len(self.x_edges))
            z_down = np.empty(len(self.x_edges))
            n_sterol = np.empty(len(self.x_edges))

            for i_x in range(bins+1):

                print(i_x)

                count_up = np.sum(ib_up == i_x)
                count_down = np.sum(ib_down == i_x)
                count_sterol = np.sum(ib_sterol == i_x)

                n_sterol[i_x] = count_sterol

                if count_up > 0:

                    zsum_up = np.sum(head_up_z[ib_up == i_x])
                    z_up[i_x] = zsum_up/count_up
                else:
                    z_up[i_x] = 0
                if count_down > 0:

                    zsum_down = np.sum(head_down_z[ib_down == i_x])
                    z_down[i_x] = zsum_down/count_down
                else:
                    z_down[i_x] = 0

            z_profiles_up.append(z_up)
            z_profiles_down.append(z_down)
            density_profiles_sterol.append(n_sterol)

        mean_z_profile_up = np.mean(z_profiles_up, axis=0)
        mean_z_profile_down = np.mean(z_profiles_down, axis=0)

        volume = 100*self.area_yz()*self.bin_width
        print(volume)

        self.mean_z_profiles_up.append(mean_z_profile_up)
        self.mean_z_profiles_down.append(mean_z_profile_down)
        mean_density_profile = (np.sum(density_profiles_sterol, axis=0)) / volume

        self.mean_density_profiles_sterol.append(mean_density_profile)

        delta_h = mean_z_profile_up - mean_z_profile_down
        self.delta_h_profiles.append(delta_h)

        y = delta_h - np.mean(delta_h)

        X1 = fft(y)

        self.fft_x.append(fftfreq(len(self.x)))
        self.fft_y.append(np.abs(X1))

    def run_fourier_transform(self,end,step):

        self.mean_density_profiles_sterol = []
        self.mean_z_profiles_up = []
        self.mean_z_profiles_down = []
        self.delta_h_profiles = []
        self.fft_x = []
        self.fft_y = []

        for i in np.arange(0,end,step):
            self.fourier_transform(i,i + step)

        a = np.array(self.delta_h_profiles)
        b = np.array(self.fft_x)
        c = np.array(self.fft_y)
        d = np.array(self.mean_z_profiles_up)
        e = np.array(self.mean_z_profiles_down)
        g = np.array(self.mean_density_profiles_sterol)

        print(d)
        print(np.shape(d))

        pp = PdfPages('z_profile_fourier.pdf')

        f = plt.figure()
        ax = f.add_subplot(111)

        for i,j in enumerate(a):

            ax.plot(self.x_edges, a[i,:], label = 'frame {}'.format(i))

        plt.legend(loc='best')
        plt.xlabel(r'X [$\AA$]')
        plt.ylabel(r'$\Delta h$ [$\AA$]')
        pp.savefig()

        f = plt.figure()
        ax = f.add_subplot(111)

        for i,j in enumerate(b):

            ax.plot(b[i,:], c[i,:], label = 'frame {}'.format(i))

        plt.legend(loc='best')
        plt.xlabel('K')
        plt.ylabel('FT')
        pp.savefig()

        f = plt.figure()
        ax = f.add_subplot(111)

        for i,j in enumerate(d):

            ax.plot(self.x_edges, d[i,:], label = 'frame {}'.format(i))
            ax.plot(self.x_edges, e[i,:], label = 'frame {}'.format(i))

        plt.legend(loc='best')
        plt.xlabel(r'X [$\AA$]')
        plt.ylabel(r'h [$\AA$]')
        pp.savefig()

        f = plt.figure()
        ax = f.add_subplot(111)
        for i,j in enumerate(g):

            ax.plot(self.x_edges, g[i,:], label = 'frame {}'.format(i))

        plt.legend(loc='best')
        plt.xlabel(r'X [$\AA$]')
        plt.ylabel('Density [number/volume]')
        pp.savefig()

        f = plt.figure()
        ax = f.add_subplot(111)

        plt.gca().set_prop_cycle(None)
        for i,j in enumerate(a):

            fit = np.polyfit(a[i,:],g[i,:],1)
            fit_fn = np.poly1d(fit)
            color = next(ax._get_lines.prop_cycler)['color']
            ax.plot(a[i,:], g[i,:], 'o', a[i,:], fit_fn(a[i,:]),'-', markeredgecolor='none', label = 'frame {}'.format(i), color=color)
            ax.plot(a[i,:], fit_fn(a[i,:]),'-', color=color)

        plt.legend(loc='best')
        plt.xlabel(r'$\Delta h$ [$\AA$]')
        plt.ylabel('Density [number/volume]')

        pp.savefig()

        pp.close()