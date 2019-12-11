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

class Order(SysInfo):


    def order_parameter_gromacs(self, start_time,end_time):

        '''This function calulate the charge density through Gromacs'''
        
        cwd = os.getcwd()
        # sn1 = ''.join(cwd + '/' + 'order_sn1.dat')
        # sn2 = ''.join(cwd + '/' + 'order_sn2.dat')

        get_selection = [GMXNAME, 'make_ndx', '-f', self.gropath, '-o', 'sn1']
        print(get_selection)
        
        inp_str = b'r DPPC & a C21\nr DPPC & a C22\nr DPPC & a C23\nr DPPC & a C24\nr DPPC & a C25\nr DPPC & a C26\nr DPPC & a C27\nr DPPC & a C28\nr DPPC & a C29\nr DPPC & a C210\nr DPPC & a C211\nr DPPC & a C212\nr DPPC & a C213\nr DPPC & a C214\nr DPPC & a C215\nr DPPC & a C216\ndel 0-4\nq\n'

        out, err = exec_gromacs(get_selection, inp_str)

        get_selection = [GMXNAME, 'make_ndx', '-f', self.gropath, '-o', 'sn2']
        
        inp_str = b'a C31\na C32\na C33\na C34\na C35\na C36\na C37\na C38\na C39\na C310\na C311\na C312\na C313\na C314\na C315\na C316\ndel 0-4\nq\n'

        out, err = exec_gromacs(get_selection, inp_str)

        get_order = [GMXNAME, 'order', '-f', self.trjpath, '-s', self.tprpath, '-n', 'sn1.ndx', '-b', str(start_time), '-e', str(end_time), \
            '-d', 'z', '-od', 'order_sn1_raw']
        
        out, err = exec_gromacs(get_order)
        
        get_order = [GMXNAME, 'order', '-f', self.trjpath, '-s', self.tprpath, '-n', 'sn2.ndx', '-b', str(start_time), '-e', str(end_time), \
            '-d', 'z', '-od', 'order_sn2_raw']
        
        out, err = exec_gromacs(get_order)
        
        with open("gmx_order.log","a") as logfile:
            logfile.write(err)
            logfile.write(out)

        with open('order_sn1_raw.xvg', 'r') as f:
            ls = f.readlines()

        with open('order_sn1.dat' , 'w') as fout:

            fout.write('Carbon_n\tScd\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
        
        with open('order_sn2_raw.xvg', 'r') as f:
            ls = f.readlines()

        with open('order_sn2.dat' , 'w') as fout:

            fout.write('Carbon_n\tScd\n')
            for l in ls:
                lc = l.strip()
                if lc:
                    if lc[0] != '#' and lc[0] != '@':
                        lf = lc.split()
                        fout.write('{}\t{}\n'.format(lf[0],lf[1]))
