'''
    This module focuses on the analysis of structural features of lipids in a bilayer

'''
import re
import os
from .. import log
from ..common import exec_gromacs, GMXNAME
from ..systeminfo import SysInfo
from ..definitions import lipidmolecules
import MDAnalysis as mda

LOGGER = log.LOGGER


def area_per_lipid(systeminfo):

    cwd = os.getcwd()
    edrout = ''.join(cwd + '/' + 'box_xy')

    box_size_arglist = [GMXNAME, 'energy', '-f', systeminfo.edrpath, '-o', edrout]
      
    inp_str = b'Box-X\nBox-Y\n0\n'
    out, err = exec_gromacs(box_size_arglist, inp_str)

    with open("gmx_box_size.log","a") as logfile:
        logfile.write(err)
        logfile.write(out)

    number_of_lipids = systeminfo.number_of_lipids/2
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


def area_per_lipid_mainlipid(systeminfo):


    u = mda.Universe(systeminfo.tprpath,systeminfo.trjpath)

    lipid_types_first = ''.join(systeminfo.molecules[0])
        
    main_lipid_selection = u.select_atoms('resname {} and name P'.format(lipid_types_first))
    number_of_main_lipid_selection = len(main_lipid_selection)
    print(number_of_main_lipid_selection)
    LOGGER.debug("number of main lipids: %d", number_of_main_lipid_selection)

    cwd = os.getcwd()
    edrout = ''.join(cwd + '/' + 'box_xy')

    box_size_arglist = [GMXNAME, 'energy', '-f', systeminfo.edrpath, '-o', edrout]
      
    inp_str = b'Box-X\nBox-Y\n0\n'
    out, err = exec_gromacs(box_size_arglist, inp_str)

    with open("gmx_box_size.log","a") as logfile:
        logfile.write(err)
        logfile.write(out)

    number_of_lipids = number_of_main_lipid_selection/2
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

